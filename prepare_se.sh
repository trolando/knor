#!/bin/bash

DIR=`dirname $0`

# build Release version of knor
cd $DIR
rm -rf build_se knor.tar.gz
cmake -B build_se -DCMAKE_BUILD_TYPE=Release .
cmake --build build_se --target knor --config Release -j 4

# make archive
cp build_se/knor se/bin
cd se
tar zcfv ../knor.tar.gz *
cd ..

# cleanup
rm -rf se/bin/knor build_se
