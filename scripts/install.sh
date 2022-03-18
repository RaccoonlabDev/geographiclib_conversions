#!/bin/bash
SCRIPT_FOLDER="$(dirname "$0")"
PACKAGE_FOLDER=$SCRIPT_FOLDER/..
WMM_FOLDER=$PACKAGE_FOLDER/wmm2020/magnetic

mkdir -p /usr/local/share/GeographicLib/magnetic                                            &&  \
    cp $WMM_FOLDER/wmm2020.wmm      /usr/local/share/GeographicLib/magnetic/wmm2020.wmm     &&  \
    cp $WMM_FOLDER/wmm2020.wmm.cof  /usr/local/share/GeographicLib/magnetic/wmm2020.wmm.cof