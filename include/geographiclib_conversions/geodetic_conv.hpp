#ifndef GEODETIC_CONVERTER_H_
#define GEODETIC_CONVERTER_H_

#include <geographiclib_conversions/geographiclib_conversions.h>

namespace geodetic_converter {

class GeodeticConverter
{
 public:

  GeodeticConverter()
  {
    haveReference_ = false;
  }

  ~GeodeticConverter()
  {
  }

  // Default copy constructor and assignment operator are OK.

  bool isInitialised()
  {
    return haveReference_;
  }

  void getReference(double* latitude, double* longitude, double* altitude)
  {
    *latitude = initial_latitude_;
    *longitude = initial_longitude_;
    *altitude = initial_altitude_;
  }

  void initialiseReference(const double latitude, const double longitude, const double altitude)
  {
    // Save NED origin
    initial_latitude_ = latitude;
    initial_longitude_ = longitude;
    initial_altitude_ = altitude;

    haveReference_ = true;
  }

  void geodetic2Enu(const double latitude, const double longitude, const double altitude,
                    double* east, double* north, double* up)
  {
    double aux_east, aux_north, aux_up;
    geographiclib_conversions::GeodeticToENU(
        latitude, longitude, altitude, aux_east, aux_north, aux_up, initial_latitude_, initial_longitude_, initial_altitude_);
    *east = aux_east;
    *north = aux_north;
    *up = aux_up;
  }

 void enu2Geodetic(const double east, const double north, const double up, 
                   double* latitude, double* longitude, double* altitude)
 {
   // Local ENU position to geodetic coordinates
   double res_latitude, res_longitude, res_altitude;
   geographiclib_conversions::ENUToGeodetic(east, north, up, 
                                            res_latitude, res_longitude, res_altitude,
                                            initial_latitude_, initial_longitude_, initial_altitude_);
    *latitude = res_latitude;
    *longitude = res_longitude;
    *altitude = res_altitude;
 }

 private:
  inline
  double rad2Deg(const double radians)
  {
    return (radians / M_PI) * 180.0;
  }

  inline
  double deg2Rad(const double degrees)
  {
    return (degrees / 180.0) * M_PI;
  }

  double initial_latitude_;
  double initial_longitude_;
  double initial_altitude_;

  bool haveReference_;

}; // class GeodeticConverter
}; // namespace geodetic_conv

#endif // GEODETIC_CONVERTER_H_
