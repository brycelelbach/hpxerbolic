#! /bin/bash

export HPX_LOCATION=~/install/hpx/gcc-4.8-release-CONTROL/
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HPX_LOCATION/lib/pkgconfig
g++-4.6 -O3 -o 1d_isothermal_flow 1d_isothermal_flow.cpp `pkg-config --cflags --libs hpx_application` -DHPX_APPLICATION_NAME=1d_isothermal_flow
#g++-4.6 -g -o 1d_isothermal_flow 1d_isothermal_flow.cpp `pkg-config --cflags --libs hpx_application_debug` -DHPX_APPLICATION_NAME=1d_isothermal_flow

#export HPX_LOCATION=~/install/hpx/gcc-4.6.2-release/
#g++-4.6 -g -o hydro_stencil hydro_stencil.cpp `pkg-config --cflags --libs hpx_application_debug` -DHPX_APPLICATION_NAME=hydro_stencil.cpp
#g++-4.6 -g -o skeleton skeleton.cpp `pkg-config --cflags --libs hpx_application_debug` -DHPX_APPLICATION_NAME=skeleton.cpp

#g++-4.6 -O3 -std=c++0x -o 1d_isothermal_flow_prehpx 1d_isothermal_flow_prehpx.cpp -I/opt/boost/1.55.0-release

#g++-4.6 -g -o 1d_isothermal_flow_hpx 1d_isothermal_flow_hpx.cpp `pkg-config --cflags --libs hpx_application_debug` -DHPX_APPLICATION_NAME=1d_isothermal_flow_hpx

#g++-4.6 -g -o 1d_isothermal_flow_cleaned 1d_isothermal_flow_cleaned.cpp `pkg-config --cflags --libs hpx_application_debug` -DHPX_APPLICATION_NAME=1d_isothermal_flow_cleaned

