#!/bin/bash

unzip bwa.zip
tar xzf Python-2.7.tgz
tar xjf samtools-1.16.1.tar.bz2

cd bwa-master
make

cd ../Python-2.7
./configure
make

cd ../samtools-1.16.1
./configure
make install