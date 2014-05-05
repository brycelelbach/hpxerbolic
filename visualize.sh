#! /bin/bash

mkdir -p scratch
cd scratch

rm step.*.jpeg > /dev/null 2>&1

SIZE=1280x800
STEPS=1000

# every INCth step is rendered
INC=1

FPS_FACTOR=25

seq 1 $INC $(($STEPS-1)) | parallel --eta 'gnuplot -e "STEP={}" -e "INPUT=\"../U_hpx.dat\"" ../plot_movie_frame.gpi'

ffmpeg -y -pattern_type glob -r $(($FPS_FACTOR/$INC)) -i './step.*.jpeg' -s $SIZE -b:v 4096k -r 24 U_hpx.mpeg
mv -f U_hpx.mpeg ..

cd ..

gnuplot plot_U_hpx.gpi
gnuplot plot_dt_hpx.gpi
