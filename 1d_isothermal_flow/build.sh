#! /bin/bash

export HPX_LOCATION=~/install/hpx/gcc-4.8-release-CONTROL/
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HPX_LOCATION/lib/pkgconfig
g++-4.6 -O3 -o 1d_isothermal_flow 1d_isothermal_flow.cpp `pkg-config --cflags --libs hpx_application` -DHPX_APPLICATION_NAME=1d_isothermal_flow
