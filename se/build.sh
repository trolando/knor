#!/bin/bash
DIR=`dirname $0`/
cd ${DIR}source
mkdir build
cd build
cmake3 ..
make
cd ../..
cp ${DIR}source/build/knor ${DIR}bin/knor
