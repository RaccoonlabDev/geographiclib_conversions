/*
 * Copyright (c) 2017, Amber Garage, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Willow Garage, Inc. nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <geographiclib_conversions/geographiclib_conversions.h>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include <GeographicLib/MagneticModel.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace GeographicLib;
using namespace std;

void geographiclib_conversions::GeodeticToENU(double latitude, double longitude, double altitude,
                                               double &east, double &north, double &up,
                                               double ref_latitude, double ref_longitude, double ref_altitude)
{
    try {
        Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
        LocalCartesian orig(ref_latitude, ref_longitude, ref_altitude, earth);

        // forward calculation
        orig.Forward(latitude, longitude, altitude,
                     east, north, up);
    }
    catch (const exception& e) {
        cerr << "Exception when converting to local coordinates in convert_global_position_to_local_position" << endl;
    }
}

void geographiclib_conversions::ENUToGeodetic(const double east, const double north, const double up, 
                                               double& latitude, double& longitude, double& altitude,
                                               const double ref_latitude, const double ref_longitude, const double ref_altitude)
{
    try {
        Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
        LocalCartesian orig(ref_latitude, ref_longitude, ref_altitude, earth);

        // forward calculation
        orig.Reverse(east, north, up,
                     latitude, longitude, altitude); 
    }
    catch (const exception& e) {
        cerr << "Exception when converting to global coordinates in convert_local_position_to_global_position" << endl;
    }
}

void geographiclib_conversions::MagneticField(const double lat, const double lon, const double h, double& Bx, double& By, double& Bz)
{
    try {
        MagneticModel mag("wmm2020");
        double t = 2020;
        double Bx_nt, By_nt, Bz_nt;
        mag(t, lat,lon, h, Bx_nt, By_nt, Bz_nt);
        // nanotesla to gauss
        Bx = Bx_nt/1e5f;
        By = By_nt/1e5f;
        Bz = Bz_nt/1e5f;
        // cout << " " << lat << " " << lon << " " << h << " " << Bx<< " "<<By<< " "<<Bz<<endl;
        // Bx the easterly component of the magnetic field (nanotesla).
        // By the northerly component of the magnetic field (nanotesla).
        // Bz the vertical (up) component of the magnetic field (nanotesla).
    }
    catch (const exception& e) {
        // sudo geographiclib-get-magnetic all
        cerr << "Caught exception: " << e.what() << "\n";
    }
}

float geographiclib_conversions::MagneticDeclination(const double latitude, const double longitude, const double altitude)
{
    try {
        MagneticModel mag("wmm2020");
        double Bx, By, Bz;
        mag(2017,
            latitude,     /*Latitude*/
            longitude,     /*Longitude*/
            altitude,     /*Altitude*/
            Bx, By, Bz);        /*[out] components of the magnetic field, enu (nanotesla)*/
        double H, F, D, I;
        MagneticModel::FieldComponents(Bx, By, Bz, H, F, D, I);
        return D;               /*the declination of the field (degrees east of north)*/
    }
    catch (const exception& e) {
        // sudo geographiclib-get-magnetic all
        cerr << "Caught exception: " << e.what() << " Did you installed WWM? Look at https://geographiclib.sourceforge.io/html/magnetic.html#magneticinst" << endl;
        return 0;
    }
}

void geographiclib_conversions::NormalizeGeodetic(double& latitude, double& longitude) 
{
    longitude = (fmod(fmod((longitude + 180.0), 360.0) + 360.0, 360.0) - 180.0);
    latitude = min(max(latitude, -90.0), 90.0);
}

bool geographiclib_conversions::IsGeodeticValid(const double latitude, const double longitude) 
{
    if (latitude < -90.0 || latitude > 90.0)
        return false;
    if (longitude < -180.0 || longitude >= 180.0)
        return false;
    return true;
}