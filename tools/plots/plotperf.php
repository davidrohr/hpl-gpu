<?php

$dgemm_peaks = array();
if (file_exists('../dgemm_peak.list'))
{
	$dgemm_tmp = explode("\n", str_replace("\r", '', file_get_contents('../dgemm_peak.list')));
	foreach ($dgemm_tmp as $line)
	{
		if (strlen($line) == 0) continue;
		$tmp = explode(" ", $line);
		$dgemm_peaks[$tmp[0]] = $tmp[1];
	}
}

$names = array();
if (file_exists('../names.txt'))
{
	$names_tmp = explode("\n", str_replace("\r", '', file_get_contents('../names.txt')));
	foreach ($names_tmp as $line)
	{
		if (strlen($line) == 0) continue;
		$tmp = explode(" ", $line, 2);
		$names[$tmp[0]] = $tmp[1];
	}
}

$nbs = array();
$effs = array();

$maxn = 0;

$result = array();

$la2_timing = 1;

$plots = array();
$power_avail = 1;
for ($i = 1;$i < count($argv);$i++)
{
	$file = file_get_contents($argv[$i] . '.filter');
	$lines = explode("\n", str_replace("\r", '', $file));
	$iteration = "0";
	$totaltime = 0;
	$totaliteration = false;
	$gflops = 0;
	$vals = array();
	$extracputime = 0;
	$extracputimeuse = 0;
	$gputime = 0;
	$cputime = 0;
	$facttime = 0;
	$dtrmtime = 0;
	$laswptime = 0;
	$bcasttime = 0;
	$ubcasttime = 0;
	$dtrsmpre = 0;
	$laswppre = 0;
	$ubcastpre = 0;
	$ubcastcount = $la2_timing ? 0 : 3;
	$dlattime = 0;
	$dlatime = 0;
	$dtrsmtime = 0;
	$usenb = 0;
	$lookahead2b = 0;
	$inwarmup = 0;
	foreach ($lines as $line)
	{
		if (strpos($line, "Running warmup iteration") === 0) $inwarmup = 1;
		else if ($inwarmup)
		{
			if (strpos($line, "Timer ITERATION") === 0) $inwarmup = 0;
			continue;
		}

		if (strpos($line, "Iteration") === 0)
		{
			$tmp = explode(" ", $line);
			$tmp = explode("=", $tmp[3]);
			$iteration = $tmp[1];
			if ($totaliteration === false) $totaliteration = $iteration;
			if ($iteration > $maxn) $maxn = $iteration;
		}
		else if (strpos($line, "Timer ITERATION") === 0)
		{
			$tmp = explode(" ", $line);
			$time = (float) $tmp[8];
			$totaltime += $time;
			$vals[] = array('n' => $iteration, 'perf' => $gflops, 'lperf' => round(min($gflops, $gflops * $fulldgemmtime / $time), 2), 'time' => $time, 'cputime' => round($cputime, 2), 'extracputime' => round($extracputimeuse, 2), 'gputime' => round($gputime, 2), 'totaltime' => $totaltime, 'nb' => $usenb, 'facttime' => round($facttime, 2), 'dtrsmtime' => round($dtrsmtime, 2), 'laswptime' => round($laswptime, 2), 'dlattime' => round($dlattime, 2), 'bcasttime' => round($bcasttime, 2), 'ubcasttime' => round($ubcasttime, 2), 'dgemmtime' => round($dgemmtime, 2), 'dtrsmpre' => round($dtrsmpre, 2), 'laswppre' => round($laswppre, 2), 'ubcastpre' => round($ubcastpre, 2), 'dlatime' => round($dlatime, 2));
			$facttime = 0;
			$dtrmtime = 0;
			$laswptime = 0;
			$bcasttime = 0;
			$ubcasttime = 0;
			$dlattime = 0;
			$dlatime = 0;
			$dtrsmtime = 0;
			$dtrsmpre = 0;
			$laswppre = 0;
			$ubcastpre = 0;
			$ubcastcount = $la2_timing ? 0 : 3;
		}
		else if (strpos($line, 'System Gflops') !== false)
		{
			$tmp = explode(" ", $line);
			$gflops = round($tmp[count($tmp) - 1], 2);
			$fulldgemmtime = (float) $tmp[count($tmp) - 4];
			$dgemmtime = $fulldgemmtime - ($dtrsmpre + $laswppre + $ubcastpre);
			$extracputime = 0;
		}
		else if (strpos($line, 'Timer RPFACT') === 0)
		{
			$tmp = explode(" ", $line);
			$extracputime += (float) $tmp[count($tmp) - 1];
			$facttime = (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, 'Timer DTRSM') === 0)
		{
			$tmp = explode(" ", $line);
			if ($dtrsmpre == 0 && $la2_timing)
				$dtrsmpre = (float) $tmp[count($tmp) - 1];
			else
				$dtrsmtime += (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, 'Timer LASWP') === 0)
		{
			$tmp = explode(" ", $line);
			if ($laswppre == 0 && $la2_timing)
				$laswppre = (float) $tmp[count($tmp) - 1];
			else
				$laswptime += (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, 'Timer BCAST') === 0)
		{
			$tmp = explode(" ", $line);
			$extracputime += (float) $tmp[count($tmp) - 1];
			$bcasttime = (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, 'Timer UBCAST') === 0)
		{
			$tmp = explode(" ", $line);
			if ($ubcastcount++ < 3)
				$ubcastpre += (float) $tmp[count($tmp) - 1];
			else
				$ubcasttime += (float) $tmp[count($tmp) - 1];
			if ($ubcastcount >= 4) $lookahead2b = 1;
		}
		else if (strpos($line, 'Timer DLATCPY') === 0)
		{
			$tmp = explode(" ", $line);
			$dlattime = (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, 'Timer DLACPY') === 0)
		{
			$tmp = explode(" ", $line);
			$dlatime = (float) $tmp[count($tmp) - 1];
			$extracputime += (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, "Starting DGEMM Run") === 0)
		{
			$extracputimeuse = $extracputime;
			$facttime = 0;
			$dtrsmtime = 0;
			$bcasttime = 0;
		}
		else if (strpos($line, "GPU Ratio") !== false)
		{
			$tmp = $line;
			while (strpos($tmp, '  ') !== false) $tmp = str_replace('  ', ' ', $tmp);
			$tmp = explode(" ", $tmp);
			$gputime = (float) $tmp[5] - ($dtrsmpre + $laswppre + $ubcastpre);
			$cputime = (float) (strpos($line, "Total CPU Time") !== false ? ($tmp[22] - $ubcastpre - $laswppre) : $tmp[10]);
		}
		else if (strpos($line, "NB     :   ") === 0)
		{
			$tmp = explode(" ", trim($line));
			$usenb = (float) $tmp[count($tmp) - 1];
		}
		else if (strpos($line, "Timer DGEMM") !== false && strpos($line, "Timer UPDATE") !== false && strpos($line, 'WC') === 0)
		{
			echo "Ignoring Unexpected String: $line\n";
		}
	}
	foreach ($vals as &$ref_val)
	{
		$ref_val['progress'] = round(100. * $ref_val['totaltime'] / $totaltime, 2);
		$ref_val['progress_size'] = round(100 * ($totaliteration - $ref_val['n']) / $totaliteration, 2);
	}

	$pname = str_replace(array(".mout", ".out", ".power", ".temp"), "", $argv[$i]) . '.pfilter';
	$eff = 0;
	if (file_exists($pname))
	{
		$filepower = file_get_contents($pname);
		$powers = explode("\n", $filepower);
		foreach ($powers as &$ref_power) $ref_power = explode(" ", $ref_power);
		foreach ($vals as &$ref_val)
		{
			$refpower = 0;
			$refpowerstart = 99999999;
			$refpowerend = -200;
			$numvals = 0;
			foreach ($powers as $power)
			{
				if (count($power) >= 3 && $power[0] >= $ref_val['totaltime'] - $ref_val['time'] - 1 && $power[0] <= $ref_val['totaltime'] + 1)
				{
					if ($power[0] - $power[2] / 2 < $refpowerstart) $refpowerstart = $power[0] - $power[2] / 2;
					if ($power[0] + $power[2] / 2 > $refpowerend) $refpowerend = $power[0] + $power[2] / 2;
					$refpower += $power[1];
					$numvals++;
				}
				if ($power[0] == 'EFF') $eff = $power[1];
			}
			if ($refpowerend != $refpowerstart) $powerval = $refpower / ($refpowerend - $refpowerstart) * $ref_val['time'];
			else $powerval = 0;

			//echo "N: $ref_val[n]   Time: $ref_val[time]   Power: $refpower   Scaled Power: $powerval   Powertime: " . ($refpowerend - $refpowerstart) . " ($refpowerend - $refpowerstart)   Powercount: $numvals\n";

			$ref_val['power'] = $powerval;
		}
	}
	else
	{
		foreach ($vals as &$ref_val)
		{
			$ref_val['power'] = 0;
		}
		$power_avail = 0;
	}

	$plots[] = $vals;
	if ($usenb == 0)
	{
		echo "WARNING: NB not found in $argv[$i]\n";
		$usenb = 2048;
	}
	$nbs[] = $usenb;
	$effs[] = $eff;

	$res = file_get_contents($argv[$i] . '.filter');
	$lines = explode("\n", str_replace("\r", '', $res));
	$result[$i] = '';
	foreach ($lines as $line)
	{
		if (strpos($line, 'WC') === 0 || strpos($line, 'WR') === 0)
		{
			if ($line{4} != 'L' || $line{6} != 'C') continue;
			while (strpos($line, '  ') !== false) $line = str_replace('  ', ' ', $line);
			$line = explode(' ', $line);
			$result[$i] = ' (' . round((float) $line[7], 1) . ' GFlop/s)';
			break;
		}
	}
}

$line = 0;
$outputpre = '';
$output100 = '';
$output0 = '';
for ($i = 1;$i < count($argv);$i++)
{
	$outputpre .= "\"$argv[$i]_n\"\t";		//1
	$outputpre .= "\"$argv[$i]_progress\"\t";	//2
	$outputpre .= "\"$argv[$i]_performance\"\t";	//3
	$outputpre .= "\"$argv[$i]_lperformance\"\t";	//4
	$outputpre .= "\"$argv[$i]_time\"\t";		//5
	$outputpre .= "\"$argv[$i]_time_cpu\"\t";	//6
	$outputpre .= "\"$argv[$i]_time_gpu\"\t";	//7
	$outputpre .= "\"$argv[$i]_serial_time_cpu\"\t";//8
	$outputpre .= "\"$argv[$i]_time_fact\"\t";	//9
	$outputpre .= "\"$argv[$i]_time_dtrsm\"\t";	//10
	$outputpre .= "\"$argv[$i]_time_laswp\"\t";	//11
	$outputpre .= "\"$argv[$i]_time_dlat\"\t";	//12
	$outputpre .= "\"$argv[$i]_time_bcast\"\t";	//13
	$outputpre .= "\"$argv[$i]_time_caldgemm\"\t";	//14
	$outputpre .= "\"$argv[$i]_time_ubcast\"\t";	//15
	$outputpre .= "\"$argv[$i]_time_dtrsmpre\"\t";	//16
	$outputpre .= "\"$argv[$i]_time_laswppre\"\t";	//17
	$outputpre .= "\"$argv[$i]_time_ubcastpre\"\t";	//18
	$outputpre .= "\"$argv[$i]_time_dla\"\t";	//19
	$outputpre .= "\"$argv[$i]_power\"\t";	//20
	$output100 .= "0\t100\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
	$output0 .= "$maxn\t" . $plots[$i - 1][0]['progress'] . "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
}
$outputpre .= "\n";
$output100 .= "\n";
$output0 .= "\n";
$output = '';
do
{
	$entryfound = 0;
	$tmpoutput = '';
	foreach ($plots as $plotnum => $plot)
	{
		if (isset($plot[$line]))
		{
			$entryfound = 1;
			$plot = $plot[$line];
			$tmpoutput .= "$plot[n]\t$plot[progress]\t$plot[perf]\t$plot[lperf]\t$plot[time]\t$plot[cputime]\t$plot[gputime]\t$plot[extracputime]\t$plot[facttime]\t$plot[dtrsmtime]\t$plot[laswptime]\t$plot[dlattime]\t$plot[bcasttime]\t$plot[dgemmtime]\t$plot[ubcasttime]\t$plot[dtrsmpre]\t$plot[laswppre]\t$plot[ubcastpre]\t$plot[dlatime]\t$plot[power]\t";
		}
		else
		{
			$tmpoutput .= "0\t100\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
		}
	}
	$line++;
	if ($entryfound) $output .= $tmpoutput . "\n";
} while ($entryfound);

$num_cols = 20;

file_put_contents("plot.dat", $outputpre . $output);

file_put_contents("plot2.dat", $outputpre . $output0 . $output . $output100);

$maxn += 2000;

$script = '';
$script2 = '';
for ($axis = 0;$axis < 4;$axis++)
{
	$xaxis = $axis & 1;
	$yscale = $axis & 2;
	$nameappend = $xaxis ? "_n" : "";
	if ($yscale) $nameappend .= "_normalized";
	
	for ($s = 0;$s < 4;$s++)
	{
		if ($s == 3 && $power_avail == 0) continue;
		if ($s < 3 && $axis >= 2) continue;
		$script .= '
		load "../pre.gnuplot"
		set style data lines
		set terminal pdf color enhanced font ' . (count($argv) < 8 ? 'font_normal' : (count($argv) < 12 ? 'font_small' : (count($argv) < 16 ? 'font_smaller' : 'font_smallest'))) . ' size 9.5in,5.5in dashed
		' . ($s == 0 ? 
			"set ylabel \"DGEMM Performance [GFlop/s]\"\nset output \"gp_perf$nameappend.pdf\"" : 
			($s == 1 ? "set ylabel \"Iteration Time [s]\"\nset output \"gp_time$nameappend.pdf\"" : 
				($s == 2 ? "set ylabel \"Time [s]\"\nset output \"gp_time2$nameappend.pdf\"" :
				($yscale ? "set ylabel \"Power Efficiency [MFlop/J]\"\nset output \"gp_power$nameappend.pdf\"" : "set ylabel \"Total Power [J]\"\nset output \"gp_power$nameappend.pdf\""))
			)
		) . '
		set xlabel "' . ($xaxis ? 'Iteration (Remaining Matrix Size)' : 'Progress [%]') . '"
		' . ($xaxis ? "set mxtics 5\nset grid ytics xtics mxtics\nset xrange [$maxn:-2000]" : "set xrange [0:102]") . '
		' . ($yscale ? "set mytics 5\nset grid mytics\n" : '') . '
		' . ($s == 0 ? 'set key bottom left Left reverse' : 'set key top right') . '
		plot';
		$num = 1;
		$dgemmpeakadded = 0;
		for ($i = 1;$i < count($argv);$i++)
		{
			$async_dlatcpy = 1;
			if ($s == 0)
			{
				if ($i == 1)
				{
					$script .= " '../empty' using 1:2 ls 1 lc rgb \"#777777\" title \"DGEMM Performance\",";
					$script .= " '../empty' using 1:2 ls 2 lc rgb \"#777777\" title \"HPL Performance\",";
				}
				foreach ($dgemm_peaks as $index => $val)
				{
					if (strpos($argv[$i], $index) === 0)
					{
						$script .= " $val ls 3 lc rgb user_color$num notitle,";
						if ($dgemmpeakadded == 0)
						{
							$script .= "'../empty' using 1:2 ls 3 lc rgb \"#777777\" title \"DGEMM Peak Performance\",";
						}
						$dgemmpeakadded = 1;
						break;
					}
				}
				$script .= "'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":" . (($i - 1)*$num_cols+4) . " ls 2 lc rgb user_color$num notitle,";
				$script .= "'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":" . (($i - 1)*$num_cols+3) . " ls 1 lc rgb user_color$num title \"" . (isset($names[$argv[$i]]) ? $names[$argv[$i]] : str_replace(array("_", '.out'), array("-"), $argv[$i])) . $result[$i] . "\"";
				if (1 || $xaxis == 0)
				{
					$scriptplot = "";
					$timesum = "0";
					$timemin = "0";
					$scriptplots = array( //col-offset, color, title, new time min
						array(15, 4, "UBCAST Time", false),
						array(11, 8, "LASWP Time", false),
						array(10, 3, "DTRSM Time", false),
						array(19, 5, "DLACPY Time", false),
						$async_dlatcpy ? array(12, 6, "DLATCPY", false) : false,
						array(9, 9, "Factorization Time", false),
						array(13, 2, "BCAST Time", false),
						array(6, 1, false, 0),	//CPU DGEMM
						array(17, 8, false, 14),
						array(18, 4, false, false),
						array(16, 3, false, false),
						!$async_dlatcpy ? array(12, 6, "DLATCPY Time", false) : false,

						array(8, 7, "Lookahead Initialization Time", false),
					);
					foreach ($scriptplots as $tmp)
					{
						if ($tmp === false) continue;
						if ($tmp[3] !== false)
						{
							if ($tmp[3] == 0)
							{
								$timesum = $timemin = "0";
							}
							else
							{
								$timesum = $timemin = '$' . (($i - 1)*$num_cols+$tmp[3]);
							}
						}
						$timesum .= ' + $' . (($i - 1)*$num_cols+$tmp[0]);
						$scriptplot = "'plot2.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":($timesum):($timemin) w filledcurves ls $tmp[1] lt rgb user_bgcolor$tmp[1] " . ($tmp[2] !== false ? ("title \"$tmp[2]\"") : "notitle") . ", \\\n$scriptplot";
					}

					$script2 .= "

						load \"../pre.gnuplot\"
						set style data lines
						set terminal pdf color enhanced font font_normal size 9.5in,5.5in solid
						set ylabel \"Time [s]\"
						
						set output \"time_$argv[$i]" . ($xaxis ? "_n" : '') . ".pdf\"
						set xlabel \"" . ($xaxis ? 'Iteration (Remaining Matrix Size)' : 'Progress [%]') . "\"
						" . ($xaxis ? "set mxtics 5\nset grid ytics xtics mxtics\nset xrange [$maxn:-2000]" : "set xrange [0:102]") . "
						#set yrange [0:25]
						set key top right invert
						plot \\
						$scriptplot 'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":" . (($i - 1)*$num_cols+6) . " lw 5 lt rgb user_color1 title \"CALDGEMM CPU Time\", \\
						'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":" . (($i - 1)*$num_cols+7) . " lw 5 lt rgb user_color10 title \"CALDGEMM GPU Time\"
					";
				}
			}
			else if ($s == 1)
			{
				$script .= "'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":(\$" . (($i - 1)*$num_cols+5) . '*' . (2048/$nbs[$i - 1]) . ") ls 1 lc rgb user_color$num title \"" . (isset($names[$argv[$i]]) ? $names[$argv[$i]] : str_replace(array("_", '.out'), array("-"), $argv[$i])) . $result[$i] . "\"";
			}
			else if ($s == 2)
			{
				if ($dgemmpeakadded == 0)
				{
					$dgemmpeakadded = 1;
					$script .= "'../empty' using 1:2 ls 1 lc rgb \"#777777\" title \"CPU Time\",";
					$script .= "'../empty' using 1:2 ls 2 lc rgb \"#777777\" title \"GPU Time\",";
				}
				$script .= "'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":(\$" . (($i - 1)*$num_cols+6) . '*' . (2048./$nbs[$i - 1]) . ") ls 1 lc rgb user_color$num title \"" . (isset($names[$argv[$i]]) ? $names[$argv[$i]] : str_replace(array("_", '.out'), array("-"), $argv[$i])) . $result[$i] . "\",";
				$script .= "'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":(\$" . (($i - 1)*$num_cols+7) . '*' . (2048./$nbs[$i - 1]) . ") ls 2 lc rgb user_color$num notitle";
			}
			else if ($s == 3)
			{
				$script .= "'plot.dat' using " . (($i - 1)*$num_cols+2-$xaxis) . ":(" . ($yscale ? (
				  "\$" . (($i - 1)*$num_cols+1) . " > 5000 ? ((1e-9 * \$" . (($i - 1)*$num_cols+1) . " * \$" . (($i - 1)*$num_cols+1) . " * (2./3. * \$" . (($i - 1)*$num_cols+1) . "+ .3/2.) - 1e-9 * (\$" . (($i - 1)*$num_cols+1) . " - " . $nbs[$i - 1] . ") * (\$" . (($i - 1)*$num_cols+1) . " - " . $nbs[$i - 1] . ") * (2./3. * (\$" . (($i - 1)*$num_cols+1) . "-" . $nbs[$i - 1] . ")+ 3./2.)) / \$" . (($i - 1)*$num_cols+20) . ") : 0"
				  ) : "\$" . (($i - 1)*$num_cols+20)) . ") ls 1 lc rgb user_color$num title \"" . (isset($names[$argv[$i]]) ? $names[$argv[$i]] : str_replace(array("_", '.out'), array("-"), $argv[$i])) . $result[$i] . " (" . $effs[$i - 1] . " MFlop/W)" . "\"";
			}
			if ($i < count($argv) - 1) $script .= ',';
			$num = $num % 12 + 1;
		}
	}
}
$script .= "\n" . $script2;

file_put_contents("make.gnuplot", $script);

?>