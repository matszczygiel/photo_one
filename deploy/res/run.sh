#!/bin/bash

rm log.out

PPATH="/home/mateusz/workspace/photo_h/build/Release/"


for e in 26.704  27.283   27.848   27.929   28.643   29.426   30.276   31.194   32.181   33.235   34.358   35.548   36.807   37.728   38.133 40.8 44.0 52.4 61.9
do
	echo $e
	sed -i -e "s/^PHOTON_EN.*/PHOTON_EN $e/g" test.inp
    
    sed -i -e "s/^GAUGE.*/GAUGE            dipole/g" test.inp
    sed -i -e "s/^FILE_OUT.*/FILE_OUT         res_dip.out/g" test.inp
	$PPATH./photo test.inp $1 >> log.out
    sed -i -e "s/^GAUGE.*/GAUGE            velocity/g" test.inp
    sed -i -e "s/^FILE_OUT.*/FILE_OUT         res_vel.out/g" test.inp
	$PPATH./photo test.inp $1 >> log.out

done