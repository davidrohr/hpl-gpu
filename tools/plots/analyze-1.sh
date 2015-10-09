#!/bin/bash
rm "`pwd`-1/"*
for i in `ls *.mout`; do
  cp `echo $i | sed "s/\.mout/-0.out/g"` "`pwd`-1/`echo $i | sed "s/\.mout/.out/g"`"
  cp `echo $i | sed "s/\.mout/.power/g"` "`pwd`-1/"
done
pushd . > /dev/null
cd "`pwd`-1"
../run.sh *.out
popd > /dev/null
