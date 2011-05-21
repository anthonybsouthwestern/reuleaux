#!/usr/bin/bash -x

VERSION="$1"
SUBVERSION="$2"

DIRNAME=reuleaux-$VERSION\.$SUBVERSION

mkdir rel/$DIRNAME
mkdir rel/$DIRNAME/{src,bin,lib,data}

cp src/*.cpp rel/$DIRNAME/src
cp src/*.h rel/$DIRNAME/src

cp src/makefile rel/$DIRNAME/src/makefile

cp makefile rel/$DIRNAME
cp README.txt rel/$DIRNAME
cp LICENSE.txt rel/$DIRNAME

DATADIRS="ReuleauxExample"

for i in $DATADIRS; do
    mkdir rel/$DIRNAME/data/$i
    cp data/$i/$i.psc rel/$DIRNAME/data/$i
    cp data/$i/$i.cfg.rel rel/$DIRNAME/data/$i/$i.cfg
done

dos2unix rel/$DIRNAME/*
dos2unix rel/$DIRNAME/src/*
for i in $DATADIRS; do
    dos2unix rel/$DIRNAME/data/$i/*
done

cd rel

tar -czf $DIRNAME.tgz $DIRNAME
