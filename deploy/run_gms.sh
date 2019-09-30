#!/bin/bash

GAMESS_PATH="/home/mateusz/gamessSTO"
PHOTO_MISC_PATH="/home/mateusz/workspace/photo_misc/build/Release/"


HERE=$(dirname $(readlink -f $0))
mkdir tmp_data
mkdir data
TMP="$HERE/tmp_data"
INPS="$HERE/inputs"
DATA="$HERE/data"

$PHOTO_MISC_PATH./main_h -gms $INPS/basis.inp $INPS

cp $INPS/h.inp $GAMESS_PATH

cd $GAMESS_PATH
./rungms h.inp 01 > log_i.out
rm h.inp

mv c.dat H_c.dat
mv e.dat H_e.dat

mv H_c.dat $DATA
mv H_e.dat $DATA
mv log_i.out $TMP

cd $TMP
grep 'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS' log_i.out > basis_l.dat
cd ..

python3 prep/format_length.py
NUM=$(cat basis_l.dat)
mv basis_l.dat basis.dat
mv basis.dat $DATA

TARGET="$HERE/res/test.inp"

sed -i -e "s/^NUMBER_GTO .*/NUMBER_GTO       ${NUM}/g" $TARGET

sed -i -e "s|^PATH_IN .*|PATH_IN       ${HERE}/|g" $TARGET

mkdir ${HERE}/res/outps
sed -i -e "s|^PATH_OUT .*|PATH_OUT       ${HERE}/res/outps/|g" $TARGET
