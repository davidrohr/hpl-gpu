#!/bin/bash
if [ -e Make.Generic ]; then
    echo Make.Generic exists, good
else
    cp setup/Make.Generic .
fi
if [ -e Make.Generic.Options ]; then
    echo Make.Generic.Options exists, good
else
    cp setup/Make.Generic.Options .
fi
make arch=Generic -j 40 $1 $2 $3 $4 $5 $6 $7 $8 $9
