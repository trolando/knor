#!/bin/bash
DIR=`dirname $0`/se
cd ${DIR}
rm -rf source ../knor.tar.gz
mkdir source
cp -r ../libs ../src ../CMakeLists.txt source
tar zcfv ../knor.tar.gz *
rm -rf source
