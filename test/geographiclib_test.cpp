#include <cmath>
#include <cstdio>

namespace AmberGarage {
namespace Trajen {
namespace GeodeticConverter {

// Geodetic system parameters
const static double kSemimajorAxis = 6378137;
const static double kSemiminorAxis = 6356752.3142;
const static double kFirstEccentricitySquared = 6.69437999014 * 0.001;
const static double kSecondEccentricitySquared = 6.73949674228 * 0.001;
const static double kFlattening = 1 / 298.257223563;


inline static
bool isValidGeodetic(const double lat, const double lng)
{
    return -90 <= lat && lat <= 90 && -180 <= lng && lng <= 180;
}


inline static
double rad2Deg(const double radians)
{
    return (radians / M_PI) * 180.0;
}

inline static
double deg2Rad(const double degrees)
{
    return (degrees / 180.0) * M_PI;
}

static 
void geodetic_to_ecef(const double latitude, const double longitude, const double altitude, 
                      double* x, double* y, double* z)
{
    // Convert geodetic coordinates to ECEF.
    // http://code.google.com/p/pysatel/source/browse/trunk/coord.py?r=22
    double lat_rad = deg2Rad(latitude);
    double lon_rad = deg2Rad(longitude);
    double xi = sqrt(1 - kFirstEccentricitySquared * sin(lat_rad) * sin(lat_rad));
    *x = (kSemimajorAxis / xi + altitude) * cos(lat_rad) * cos(lon_rad);
    *y = (kSemimajorAxis / xi + altitude) * cos(lat_rad) * sin(lon_rad);
    *z = (kSemimajorAxis / xi * (1 - kFirstEccentricitySquared) + altitude) * sin(lat_rad);
}

static 
void ecef_to_geodetic(const double x, const double y, const double z, 
                             double* latitude, double* longitude, double* altitude)
{
    // Convert ECEF coordinates to geodetic coordinates.
    // J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates
    // to geodetic coordinates," IEEE Transactions on Aerospace and
    // Electronic Systems, vol. 30, pp. 957-961, 1994.

    double r = sqrt(x * x + y * y);
    double Esq = kSemimajorAxis * kSemimajorAxis - kSemiminorAxis * kSemiminorAxis;
    double F = 54 * kSemiminorAxis * kSemiminorAxis * z * z;
    double G = r * r + (1 - kFirstEccentricitySquared) * z * z - kFirstEccentricitySquared * Esq;
    double C = (kFirstEccentricitySquared * kFirstEccentricitySquared * F * r * r) / pow(G, 3);
    double S = cbrt(1 + C + sqrt(C * C + 2 * C));
    double P = F / (3 * pow((S + 1 / S + 1), 2) * G * G);
    double Q = sqrt(1 + 2 * kFirstEccentricitySquared * kFirstEccentricitySquared * P);
    double r_0 = -(P * kFirstEccentricitySquared * r) / (1 + Q)
            + sqrt(
                    0.5 * kSemimajorAxis * kSemimajorAxis * (1 + 1.0 / Q)
                            - P * (1 - kFirstEccentricitySquared) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r);
    double U = sqrt(pow((r - kFirstEccentricitySquared * r_0), 2) + z * z);
    double V = sqrt(
            pow((r - kFirstEccentricitySquared * r_0), 2) + (1 - kFirstEccentricitySquared) * z * z);
    double Z_0 = kSemiminorAxis * kSemiminorAxis * z / (kSemimajorAxis * V);
    *altitude = U * (1 - kSemiminorAxis * kSemiminorAxis / (kSemimajorAxis * V));
    *latitude = rad2Deg(atan((z + kSecondEccentricitySquared * Z_0) / r));
    *longitude = rad2Deg(atan2(y, x));
}

static 
void ecef_to_ned(const double x, const double y, const double z, 
                 const double ref_latitude, const double ref_longitude, const double ref_altitude,
                 double* north, double* east, double* down)
{
    // Converts ECEF coordinate position into local-tangent-plane NED.
    // Coordinates relative to given ECEF coordinate frame.
    double ref_latitude_rad = deg2Rad(ref_latitude); 
    double ref_longitude_rad = deg2Rad(ref_longitude); 
    double ref_ecef_x, ref_ecef_y, ref_ecef_z;
    geodetic_to_ecef(ref_latitude, ref_longitude, ref_altitude, &ref_ecef_x, &ref_ecef_y, &ref_ecef_z);
    double dx = x - ref_ecef_x, dy = y - ref_ecef_y, dz = z - ref_ecef_z;

    double phiP = atan2(ref_ecef_z, sqrt(pow(ref_ecef_x, 2) + pow(ref_ecef_y, 2)));

    double sLat = sin(phiP), cLat = cos(phiP);
    double sLon = sin(ref_longitude_rad), cLon = cos(ref_longitude_rad);
    *north = -sLat * cLon * dx - sLat * sLon * dy + cLat * dz;
    *east = -sLon * dx + cLon * dy;
    *down = -cLat * cLon * dx - cLat * sLon * dy - sLat * dz;
}

static 
void ned_to_ecef(const double north, const double east, const double down, 
                 const double ref_latitude, const double ref_longitude, const double ref_altitude,
                 double* x, double* y, double* z)
{
    // NED (north/east/down) to ECEF coordinates
    double ref_latitude_rad = deg2Rad(ref_latitude); 
    double ref_longitude_rad = deg2Rad(ref_longitude); 
    double ref_ecef_x, ref_ecef_y, ref_ecef_z; 
    geodetic_to_ecef(ref_latitude, ref_longitude, ref_altitude, &ref_ecef_x, &ref_ecef_y, &ref_ecef_z);

    double sLat = sin(ref_latitude_rad), cLat = cos(ref_latitude_rad);
    double sLon = sin(ref_longitude_rad), cLon = cos(ref_longitude_rad);

    *x = -sLat * cLon * north - sLon * east - cLat * cLon * down + ref_ecef_x;
    *y = -sLat * sLon * north + cLon * east - cLat * sLon * down + ref_ecef_y;
    *z = cLat * north - sLat * down + ref_ecef_z;
}

static 
void geodetic_to_ned(const double latitude, const double longitude, const double altitude,
                     const double ref_latitude, const double ref_longitude, const double ref_altitude,
                     double* north, double* east, double* down)
{
    // Geodetic position to local NED frame
    double x, y, z;
    geodetic_to_ecef(latitude, longitude, altitude, &x, &y, &z);
    ecef_to_ned(x, y, z, ref_latitude, ref_longitude, ref_altitude, north, east, down);
}

static 
void ned_to_geodetic(const double north, const double east, const double down, 
                     const double ref_latitude, const double ref_longitude, const double ref_altitude,
                     double* latitude, double* longitude, double* altitude)
{
    // Local NED position to geodetic coordinates
    double x, y, z;
    ned_to_ecef(north, east, down, ref_latitude, ref_longitude, ref_altitude, &x, &y, &z);
    ecef_to_geodetic(x, y, z, latitude, longitude, altitude);
}

static 
void geodetic_to_enu(const double latitude, const double longitude, const double altitude,
                     const double ref_latitude, const double ref_longitude, const double ref_altitude,
                     double* east, double* north, double* up)
{
    // Geodetic position to local ENU frame
    double x, y, z;
    geodetic_to_ecef(latitude, longitude, altitude, &x, &y, &z);

    double aux_north, aux_east, aux_down;
    ecef_to_ned(x, y, z, ref_latitude, ref_longitude, ref_altitude, &aux_north, &aux_east, &aux_down);

    *east = aux_east;
    *north = aux_north;
    *up = -aux_down;
}

static 
void enu_to_geodetic(const double east, const double north, const double up, 
                     const double ref_latitude, const double ref_longitude, const double ref_altitude,
                     double* latitude, double* longitude, double* altitude)
{
    // Local ENU position to geodetic coordinates

    const double aux_north = north;
    const double aux_east = east;
    const double aux_down = -up;
    double x, y, z;
    ned_to_ecef(aux_north, aux_east, aux_down, ref_latitude, ref_longitude, ref_altitude, &x, &y, &z);
    ecef_to_geodetic(x, y, z, latitude, longitude, altitude);
}

static
void geodetic_to_enu_simple(const double latitude, const double longitude, const double altitude,
                            const double ref_latitude, const double ref_longitude, const double ref_altitude,
                            double* east, double* north, double* up)
{
    double dx = longitude - ref_longitude;
    double dy = latitude - ref_latitude;
    double dz = altitude - ref_altitude;
    *east = deg2Rad(dx) * kSemimajorAxis * std::cos(deg2Rad(ref_latitude));
    *north = deg2Rad(dy) * kSemimajorAxis;
    *up = dz;
}

static
void enu_to_geodetic_simple(const double east, const double north, const double up,
                            const double ref_latitude, const double ref_longitude, const double ref_altitude,
                            double* latitude, double* longitude, double* altitude)
{
    *latitude = ref_latitude + rad2Deg(north / kSemimajorAxis);
    *longitude = ref_longitude + rad2Deg(east / (kSemimajorAxis * std::cos(deg2Rad(ref_latitude))));
    *altitude = ref_altitude + up;
}


// ================================================================================
// Magnetic declination
static void wmm_initial_constants(int maxord, double c[13][13], double cd[13][13], double k[13][13]);

static
void wmm_magnetic_declination(const double latitude,    // latitude in decimal degrees
                              const double longitude,   // longitude in decimal degrees
                              const double altitude,    // altitude in meters
                              const double time,        // time in decimal years, year + days/365.0
                              double *dec,              // output declination in decimal degrees
                              double *dip,              // output dip in decimal degrees
                              double *ti,               // output ti in decimal degrees
                              double *gv)               // output grid variation in decimal degrees
{
    // Initial constants
    const int MAX_ORD = 12;
    const double EPS = 1e-9;
    double c[13][13], cd[13][13], k[13][13];
    
    wmm_initial_constants(MAX_ORD, c, cd, k);

    // Start time: 2015.0. The model is valid from 2015.0 to 2020.0
    double epoch = 2015.0;
    double alt = altitude / 1000;

    double glat = latitude, glon = longitude;
    double rlat = deg2Rad(glat), rlon = deg2Rad(glon);
    double srlon = sin(rlon), srlat = sin(rlat), crlon = cos(rlon), crlat = cos(rlat);
    double srlat2 = srlat * srlat, crlat2 = crlat * crlat;
    
    /* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
    double a = 6378.137, b = 6356.7523142, re = 6371.2;
    double a2 = a * a, b2 = b * b, c2 = a2 - b2, a4 = a2 * a2, b4 = b2 * b2, c4 = a4 - b4;
    double q = sqrt(a2 - c2 * srlat2);
    double q1 = alt * q;
    double q2 = ((q1 + a2) / (q1 + b2)) * ((q1 + a2) / (q1 + b2));
    double ct = srlat / sqrt(q2 * crlat2 + srlat2);
    double st = sqrt(1.0 - (ct * ct));
    double r2 = (alt * alt) + 2.0 * q1 + (a4 - c4 * srlat2) / (q * q);
    double r = sqrt(r2);
    double d = sqrt(a2 * crlat2 + b2 * srlat2);
    double ca = (alt + d) / r;
    double sa = c2 * crlat * srlat / (r * d);
   
    double  sp[13], cp[13]; 
    sp[1] = srlon; cp[0] = 1.0; cp[1] = crlon;
    for (int m = 2; m <= MAX_ORD; m++) {
        sp[m] = sp[1] * cp[m-1] + cp[1] * sp[m-1];
        cp[m] = cp[1] * cp[m-1] - sp[1] * sp[m-1];
    }
    
    double fn[13] = {0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    double fm[13] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    double tc[13][13], pp[13], p[13][13], dp[13][13];
    pp[0] = 1.0; p[0][0] = 1.0;
    double aor = re / r;
    double ar = aor * aor;
    double br = 0.0, bt = 0.0, bp = 0.0, bpp = 0.0;
    double dt = time - epoch;
    for (int n = 1; n <= MAX_ORD; n++) {
        ar = ar * aor;
        int m = 0;
        int D4 = (n+m+1);
        while (D4 > 0) {
            /*
               COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
               AND DERIVATIVES VIA RECURSION RELATIONS
            */
            if (n == m) {
                p[m][n] = st * p[m-1][n-1];
                dp[m][n] = st * dp[m-1][n-1] + ct * p[m-1][n-1];
            } else if (n == 1 && m == 0) {
                p[m][n] = ct * p[m][n-1];
                dp[m][n] = ct * dp[m][n-1] - st * p[m][n-1];
            } else if (n > 1 && n != m) {
                if (m > n-2) {
                    p[m][n-2] = 0;
                    dp[m][n-2] = 0;
                }
                p[m][n] = ct * p[m][n-1] - k[m][n] * p[m][n-2];
                dp[m][n] = ct * dp[m][n-1] - st * p[m][n-1] - k[m][n] * dp[m][n-2];
            }

            /*
               TIME ADJUST THE GAUSS COEFFICIENTS
            */
            tc[m][n] = c[m][n] + dt * cd[m][n];
            if (m != 0) 
                tc[n][m-1] = c[n][m-1] + dt * cd[n][m-1];

            /*
               ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
            */
            double par = ar * p[m][n];
            double temp1, temp2;
            if (m == 0) {
                temp1 = tc[m][n] * cp[m];
                temp2 = tc[m][n] * sp[m];
            } else {
                temp1 = tc[m][n] * cp[m] + tc[n][m-1] * sp[m];
                temp2 = tc[m][n] * sp[m] - tc[n][m-1] * cp[m];
            }
            bt -= ar * temp1 * dp[m][n];
            bp += (fm[m] * temp2 * par);
            br += (fn[n] * temp1 * par);

            /*
               SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
            */
            if (fabs(st) < EPS && m == 1) {
                if (n == 1) 
                    pp[n] = pp[n-1];
                else 
                    pp[n] = ct * pp[n-1] - k[m][n] * pp[n-2];
                double parp = ar * pp[n];
                bpp += (fm[m] * temp2 * parp);
            }

            D4--;
            m++;
        }
    }

    if (fabs(st) < EPS)
        bp = bpp;
    else 
        bp /= st;

    /*
       ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
       GEODETIC COORDINATES
    */
    double bx = -bt * ca - br * sa;
    double by = bp;
    double bz = bt * sa - br * ca;
    //        printf("X: %.8f, Y: %.8f, Z: %.8f\n", bx, by, bz);

    /*
       COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
       TOTAL INTENSITY (TI)
    */
    double bh = sqrt((bx * bx) + (by * by));
    if (ti) *ti = sqrt((bh * bh) + (bz * bz));
    if (dec) *dec = rad2Deg(atan2(by,bx));
    if (dip) *dip = rad2Deg(atan2(bz,bh));

    /*
       COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
       GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
       (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)

       OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
    */
    if (gv) {
        *gv = -999.0;
        if (fabs(glat) >= 55.) {
            if (glat > 0.0 && glon >= 0.0) *gv = *dec - glon;
            if (glat > 0.0 && glon < 0.0) *gv = *dec + fabs(glon);
            if (glat < 0.0 && glon >= 0.0) *gv = *dec + glon;
            if (glat < 0.0 && glon < 0.0) *gv = *dec - fabs(glon);
            if (*gv > +180.0) *gv -= 360.0;
            if (*gv < -180.0) *gv += 360.0;
        }
    }
    return;

}


static
void wmm_initial_constants(int maxord, double c[13][13], double cd[13][13], double k[13][13])
{
    /*
     * The coefficients in the file follows format as:
     * 
     *      n,m,gnm,hnm,dgnm,dhnm
     * 
     * The coefficients can be loaded using code as follows:
        
     *  if (m <= n) {
     *      c[m][n] = gnm;
     *      cd[m][n] = dgnm;
     *      if (m != 0) {
     *          c[n][m-1] = hnm;
     *          cd[n][m-1] = dhnm;
     *      }
     *  }
    */
    c[0][1] = -29438.5; cd[0][1] = 10.7;
    c[1][1] = -1501.1; cd[1][1] = 17.9;
    c[1][0] = 4796.2; cd[1][0] = -26.8;
    c[0][2] = -2445.3; cd[0][2] = -8.6;
    c[1][2] = 3012.5; cd[1][2] = -3.3;
    c[2][0] = -2845.6; cd[2][0] = -27.1;
    c[2][2] = 1676.6; cd[2][2] = 2.4;
    c[2][1] = -642.0; cd[2][1] = -13.3;
    c[0][3] = 1351.1; cd[0][3] = 3.1;
    c[1][3] = -2352.3; cd[1][3] = -6.2;
    c[3][0] = -115.3; cd[3][0] = 8.4;
    c[2][3] = 1225.6; cd[2][3] = -0.4;
    c[3][1] = 245.0; cd[3][1] = -0.4;
    c[3][3] = 581.9; cd[3][3] = -10.4;
    c[3][2] = -538.3; cd[3][2] = 2.3;
    c[0][4] = 907.2; cd[0][4] = -0.4;
    c[1][4] = 813.7; cd[1][4] = 0.8;
    c[4][0] = 283.4; cd[4][0] = -0.6;
    c[2][4] = 120.3; cd[2][4] = -9.2;
    c[4][1] = -188.6; cd[4][1] = 5.3;
    c[3][4] = -335.0; cd[3][4] = 4.0;
    c[4][2] = 180.9; cd[4][2] = 3.0;
    c[4][4] = 70.3; cd[4][4] = -4.2;
    c[4][3] = -329.5; cd[4][3] = -5.3;
    c[0][5] = -232.6; cd[0][5] = -0.2;
    c[1][5] = 360.1; cd[1][5] = 0.1;
    c[5][0] = 47.4; cd[5][0] = 0.4;
    c[2][5] = 192.4; cd[2][5] = -1.4;
    c[5][1] = 196.9; cd[5][1] = 1.6;
    c[3][5] = -141.0; cd[3][5] = 0.0;
    c[5][2] = -119.4; cd[5][2] = -1.1;
    c[4][5] = -157.4; cd[4][5] = 1.3;
    c[5][3] = 16.1; cd[5][3] = 3.3;
    c[5][5] = 4.3; cd[5][5] = 3.8;
    c[5][4] = 100.1; cd[5][4] = 0.1;
    c[0][6] = 69.5; cd[0][6] = -0.5;
    c[1][6] = 67.4; cd[1][6] = -0.2;
    c[6][0] = -20.7; cd[6][0] = 0.0;
    c[2][6] = 72.8; cd[2][6] = -0.6;
    c[6][1] = 33.2; cd[6][1] = -2.2;
    c[3][6] = -129.8; cd[3][6] = 2.4;
    c[6][2] = 58.8; cd[6][2] = -0.7;
    c[4][6] = -29.0; cd[4][6] = -1.1;
    c[6][3] = -66.5; cd[6][3] = 0.1;
    c[5][6] = 13.2; cd[5][6] = 0.3;
    c[6][4] = 7.3; cd[6][4] = 1.0;
    c[6][6] = -70.9; cd[6][6] = 1.5;
    c[6][5] = 62.5; cd[6][5] = 1.3;
    c[0][7] = 81.6; cd[0][7] = 0.2;
    c[1][7] = -76.1; cd[1][7] = -0.2;
    c[7][0] = -54.1; cd[7][0] = 0.7;
    c[2][7] = -6.8; cd[2][7] = -0.4;
    c[7][1] = -19.4; cd[7][1] = 0.5;
    c[3][7] = 51.9; cd[3][7] = 1.3;
    c[7][2] = 5.6; cd[7][2] = -0.2;
    c[4][7] = 15.0; cd[4][7] = 0.2;
    c[7][3] = 24.4; cd[7][3] = -0.1;
    c[5][7] = 9.3; cd[5][7] = -0.4;
    c[7][4] = 3.3; cd[7][4] = -0.7;
    c[6][7] = -2.8; cd[6][7] = -0.9;
    c[7][5] = -27.5; cd[7][5] = 0.1;
    c[7][7] = 6.7; cd[7][7] = 0.3;
    c[7][6] = -2.3; cd[7][6] = 0.1;
    c[0][8] = 24.0; cd[0][8] = 0.0;
    c[1][8] = 8.6; cd[1][8] = 0.1;
    c[8][0] = 10.2; cd[8][0] = -0.3;
    c[2][8] = -16.9; cd[2][8] = -0.5;
    c[8][1] = -18.1; cd[8][1] = 0.3;
    c[3][8] = -3.2; cd[3][8] = 0.5;
    c[8][2] = 13.2; cd[8][2] = 0.3;
    c[4][8] = -20.6; cd[4][8] = -0.2;
    c[8][3] = -14.6; cd[8][3] = 0.6;
    c[5][8] = 13.3; cd[5][8] = 0.4;
    c[8][4] = 16.2; cd[8][4] = -0.1;
    c[6][8] = 11.7; cd[6][8] = 0.2;
    c[8][5] = 5.7; cd[8][5] = -0.2;
    c[7][8] = -16.0; cd[7][8] = -0.4;
    c[8][6] = -9.1; cd[8][6] = 0.3;
    c[8][8] = -2.0; cd[8][8] = 0.3;
    c[8][7] = 2.2; cd[8][7] = 0.0;
    c[0][9] = 5.4; cd[0][9] = 0.0;
    c[1][9] = 8.8; cd[1][9] = -0.1;
    c[9][0] = -21.6; cd[9][0] = -0.2;
    c[2][9] = 3.1; cd[2][9] = -0.1;
    c[9][1] = 10.8; cd[9][1] = -0.1;
    c[3][9] = -3.1; cd[3][9] = 0.4;
    c[9][2] = 11.7; cd[9][2] = -0.2;
    c[4][9] = 0.6; cd[4][9] = -0.5;
    c[9][3] = -6.8; cd[9][3] = 0.1;
    c[5][9] = -13.3; cd[5][9] = -0.2;
    c[9][4] = -6.9; cd[9][4] = 0.1;
    c[6][9] = -0.1; cd[6][9] = 0.1;
    c[9][5] = 7.8; cd[9][5] = 0.0;
    c[7][9] = 8.7; cd[7][9] = 0.0;
    c[9][6] = 1.0; cd[9][6] = -0.2;
    c[8][9] = -9.1; cd[8][9] = -0.2;
    c[9][7] = -3.9; cd[9][7] = 0.4;
    c[9][9] = -10.5; cd[9][9] = -0.1;
    c[9][8] = 8.5; cd[9][8] = 0.3;
    c[0][10] = -1.9; cd[0][10] = 0.0;
    c[1][10] = -6.5; cd[1][10] = 0.0;
    c[10][0] = 3.3; cd[10][0] = 0.1;
    c[2][10] = 0.2; cd[2][10] = -0.1;
    c[10][1] = -0.3; cd[10][1] = -0.1;
    c[3][10] = 0.6; cd[3][10] = 0.3;
    c[10][2] = 4.6; cd[10][2] = 0.0;
    c[4][10] = -0.6; cd[4][10] = -0.1;
    c[10][3] = 4.4; cd[10][3] = 0.0;
    c[5][10] = 1.7; cd[5][10] = -0.1;
    c[10][4] = -7.9; cd[10][4] = -0.2;
    c[6][10] = -0.7; cd[6][10] = -0.1;
    c[10][5] = -0.6; cd[10][5] = 0.1;
    c[7][10] = 2.1; cd[7][10] = 0.0;
    c[10][6] = -4.1; cd[10][6] = -0.1;
    c[8][10] = 2.3; cd[8][10] = -0.2;
    c[10][7] = -2.8; cd[10][7] = -0.2;
    c[9][10] = -1.8; cd[9][10] = -0.1;
    c[10][8] = -1.1; cd[10][8] = 0.1;
    c[10][10] = -3.6; cd[10][10] = -0.2;
    c[10][9] = -8.7; cd[10][9] = -0.1;
    c[0][11] = 3.1; cd[0][11] = 0.0;
    c[1][11] = -1.5; cd[1][11] = 0.0;
    c[11][0] = -0.1; cd[11][0] = 0.0;
    c[2][11] = -2.3; cd[2][11] = -0.1;
    c[11][1] = 2.1; cd[11][1] = 0.1;
    c[3][11] = 2.1; cd[3][11] = 0.1;
    c[11][2] = -0.7; cd[11][2] = 0.0;
    c[4][11] = -0.9; cd[4][11] = 0.0;
    c[11][3] = -1.1; cd[11][3] = 0.1;
    c[5][11] = 0.6; cd[5][11] = 0.0;
    c[11][4] = 0.7; cd[11][4] = 0.0;
    c[6][11] = -0.7; cd[6][11] = 0.0;
    c[11][5] = -0.2; cd[11][5] = 0.0;
    c[7][11] = 0.2; cd[7][11] = 0.0;
    c[11][6] = -2.1; cd[11][6] = 0.1;
    c[8][11] = 1.7; cd[8][11] = 0.0;
    c[11][7] = -1.5; cd[11][7] = 0.0;
    c[9][11] = -0.2; cd[9][11] = 0.0;
    c[11][8] = -2.5; cd[11][8] = -0.1;
    c[10][11] = 0.4; cd[10][11] = -0.1;
    c[11][9] = -2.0; cd[11][9] = 0.0;
    c[11][11] = 3.5; cd[11][11] = -0.1;
    c[11][10] = -2.3; cd[11][10] = -0.1;
    c[0][12] = -2.0; cd[0][12] = 0.1;
    c[1][12] = -0.3; cd[1][12] = 0.0;
    c[12][0] = -1.0; cd[12][0] = 0.0;
    c[2][12] = 0.4; cd[2][12] = 0.0;
    c[12][1] = 0.5; cd[12][1] = 0.0;
    c[3][12] = 1.3; cd[3][12] = 0.1;
    c[12][2] = 1.8; cd[12][2] = -0.1;
    c[4][12] = -0.9; cd[4][12] = -0.1;
    c[12][3] = -2.2; cd[12][3] = 0.0;
    c[5][12] = 0.9; cd[5][12] = 0.0;
    c[12][4] = 0.3; cd[12][4] = 0.0;
    c[6][12] = 0.1; cd[6][12] = 0.1;
    c[12][5] = 0.7; cd[12][5] = 0.0;
    c[7][12] = 0.5; cd[7][12] = 0.0;
    c[12][6] = -0.1; cd[12][6] = 0.0;
    c[8][12] = -0.4; cd[8][12] = 0.0;
    c[12][7] = 0.3; cd[12][7] = 0.0;
    c[9][12] = -0.4; cd[9][12] = 0.0;
    c[12][8] = 0.2; cd[12][8] = 0.0;
    c[10][12] = 0.2; cd[10][12] = 0.0;
    c[12][9] = -0.9; cd[12][9] = 0.0;
    c[11][12] = -0.9; cd[11][12] = 0.0;
    c[12][10] = -0.2; cd[12][10] = 0.0;
    c[12][12] = 0.0; cd[12][12] = 0.0;
    c[12][11] = 0.7; cd[12][11] = 0.0;

    double snorm[13][13];
    snorm[0][0] = 1.0, k[1][1] = 1.0;

    for (int n = 1; n <= maxord; n++) {
        snorm[0][n] = snorm[0][n-1] * (2.0*n - 1) / n;
        double j = 2.0;
        int m = 0;
        int D2 = (n - m + 1);
        while (D2 > 0) {
            k[m][n] = (((n-1) * (n-1)) - m * m) / ((2.0 * n - 1) * (2.0 * n - 3.0));
            if (m > 0) {
                double flnmj = ((n - m + 1.0) * j) / (n + m);
                snorm[m][n] = snorm[m-1][n] * sqrt(flnmj);
                j = 1.0;
                c[n][m-1] = snorm[m][n] * c[n][m-1];
                cd[n][m-1] = snorm[m][n] * cd[n][m-1];
            }
            c[m][n] = snorm[m][n] * c[m][n];
            cd[m][n] = snorm[m][n] * cd[m][n];
            D2--;
            m++;
        }
    }
}

} 
}
}

#include <geographiclib_conversions/geographiclib_conversions.h>
#include <gtest/gtest.h>
using namespace geographiclib_conversions;

TEST(GeodeticConvertTest, global2local)
{
	double latitude = 37.459749;
	double longitude = -122.209861;
	double altitude = 100;

	double ref_latitude = 37.454749;
	double ref_longitude = -122.2239861;
	double ref_altitude = 300;

	double e, n, u;
	AmberGarage::Trajen::GeodeticConverter::geodetic_to_enu(latitude, longitude, altitude, &e, &n, &u, ref_latitude, ref_longitude, ref_altitude);

	double e2, n2, u2;
	
	GeodeticToENU(latitude, longitude, altitude, e2, n2, u2, ref_latitude, ref_longitude, ref_altitude);

	EXPECT_EQ(e, e2);
	EXPECT_EQ(n, n2);
	EXPECT_EQ(u, u2);	
}

 int main(int argc, char **argv)
 {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
 }

