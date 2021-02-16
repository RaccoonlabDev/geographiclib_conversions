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

/*
 * Author: Roman Fedorenko
 */
#ifndef GEOGRAPHICLIB_CONVERSIONS_H
#define GEOGRAPHICLIB_CONVERSIONS_H

namespace geographiclib_conversions {

void GeodeticToENU(double latitude, double longitude, double altitude,
                           double& east, double& north, double& up,
                           double ref_latitude, double ref_longitude, double ref_altitude);
void ENUToGeodetic(const double east, const double north, const double up, 
                           double& latitude, double& longitude, double& altitude,
                           const double ref_latitude, const double ref_longitude, const double ref_altitude);
float MagneticDeclination(const double latitude, const double longitude, const double altitude);

void NormalizeGeodetic(double& latitude, double& longitude);

bool IsGeodeticValid(const double latitude, const double longitude);

void MagneticField(const double latitude, const double longitude, const double altitude, double& east, double& north, double& up);
}

#endif
