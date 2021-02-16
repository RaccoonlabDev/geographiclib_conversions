#include <cmath>
#include <cstdio>
#include <geographiclib_conversions/geographiclib_conversions.h>

using namespace geographiclib_conversions;

 int main(int argc, char **argv)
{
	double latitude = 37.459749;
	double longitude = -122.209861;
	double altitude = 100;

	double ref_latitude = 37.454749;
	double ref_longitude = -122.2239851;
	double ref_altitude = 300;

	double e2, n2, u2;
	
	GeodeticToENU(latitude, longitude, altitude, e2, n2, u2, ref_latitude, ref_longitude, ref_altitude);

	printf("%10.f, %10.f, %10.f\n", e2, n2, u2);

	return 0;
}