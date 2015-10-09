#!/bin/bash
rm -f plot*.dat

if [ -f /cygdrive/c/Utility/Internet/Apache/php/php.exe ]; then
    PHP=/cygdrive/c/Utility/Internet/Apache/php/php.exe
    GNUPLOT=gnuplotm.bat
else
    PHP=php
    GNUPLOT=gnuplot
fi

$PHP -f ../plotpower.php `echo $@ | grep -v filter`
$GNUPLOT make.gnuplot
RUN=
rm -f ../*.filter
for i in `echo $@ | grep -v filter`; do
	#echo $i
	MULTI=`echo $i | grep "\.mout" | wc -l`
	if [ "0$MULTI" = "01" ]; then
		PP=`cat $i | grep "\[1,0]<stdout>:P      :  " | awk '{print $3}'`
		QQ=`cat $i | grep "\[1,0]<stdout>:Q      :  " | awk '{print $3}'`
		NUMN="$(($PP * $QQ - 1))"
		for j in `seq 0 $NUMN`; do
			cat $i | grep "\[1,$j\]<\|NB     :   \|WC..L.C..   " | sed "s/\[1,[0-9]*]<std...>://" > `echo $i | sed "s/\.mout/-$j.out/"`
			RUN+=`echo $i | sed "s/\.mout/-$j.out/"`
			RUN+=" "
			cat `echo $i | sed "s/\.mout/-$j.out/"` | grep -i "Iteration\|System Gflops\|Timer \|GPU Ratio\|Starting DGEMM Run\|NB     :   \|WC..L.C..   " > `echo $i | sed "s/\.mout/-$j.out/"`.filter;
		done
	else
		cat $i | grep -i "Iteration\|System Gflops\|Timer \|GPU Ratio\|Starting DGEMM Run\|NB     :   \|WC..L.C..   " > $i.filter;
		RUN+="$i "
	fi
done

$PHP -f ../plotperf.php $RUN
$GNUPLOT make.gnuplot
rm -f *.filter *.pfilter plot*.dat make.gnuplot
