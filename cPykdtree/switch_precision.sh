#!/bin/bash

sed -i "s/$1/$2/g" src/kdtree.c
sed -i "s/$1/$2/g" src/kdtree.h
sed -i "s/$1/$2/g" src/argsort.c
sed -i "s/$1/$2/g" src/argsort.h


