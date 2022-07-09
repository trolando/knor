#!/bin/bash
DIR=`dirname $0`/
pushd ${DIR}source
cmake3 -B build -DCMAKE_BUILD_TYPE=Release .
cmake3 --build build --target knor --config Release -j 4
popd
cp ${DIR}source/build/knor ${DIR}bin/knor
