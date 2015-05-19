<?php

date_default_timezone_set('Europe/Berlin');

$script = '';

for ($j = 1;$j < count($argv);$j++)
{
	if (file_exists(str_replace(array(".mout", ".out", ".power", ".temp"), ".power", $argv[$j])))
	{
		$result = 0;
		$processingtime = -1;
		$numnodes = 1;
		$measurednodes = 1;
		if (isset($forcestarttime)) unset($forcestarttime);
		if (isset($forceendtime)) unset($forceendtime);
		if (file_exists($argv[$j]))
		{
			$res = file_get_contents($argv[$j]);
			$lines = explode("\n", str_replace("\r", '', $res));
			foreach ($lines as $line)
			{
				if (strpos($line, 'WC') !== false)
				{
					$line = substr($line, strpos($line, 'WC'));
					if ($line{4} != 'L' || $line{6} != 'C') continue;
					$line = str_replace('[1,0]<stdout>:', '', $line);
					while (strpos($line, '  ') !== false) $line = str_replace('  ', ' ', $line);
					$line = explode(' ', $line);
					$result = (float) $line[7];
					$processingtime = $line[5];
					$numnodes = $line[3] * $line[4];
				}
				else if (strpos($line, 'Calculation Start Timestamp:') !== false)
				{
					$line = str_replace('[1,0]<stdout>:', '', $line);
					$line = explode(' ', $line);
					$forcestarttime = $line[3] + 2;
				}
				else if (strpos($line, 'Calculation End Timestamp:') !== false)
				{
					$line = str_replace('[1,0]<stdout>:', '', $line);
					$line = explode(' ', $line);
					$forceendtime = $line[3];
				}
			}
		}
		//if ($numnodes) $measurednodes = 2;

		$plot = array();
		$file = file_get_contents(str_replace(array(".mout", ".out"), ".power", $argv[$j]));
		$powerout = array();
		$lines = explode("\n", str_replace("\r", '', $file));

		$min = 10000;
		$max = 0;
		$starttime = 0;
		for ($i = 1;$i < count($lines);$i++)
		{
			$line = explode(" ", $lines[$i]);
			if (count($line) < 12) continue;
			$power = (float) $line[12];
			if ($power < $min) $min = $power;
			if ($power > $max) $max = $power;
			if ($starttime == 0) $starttime = $line[0];
			$endtime = $line[0];
		}

		$mid = ($max + $min) / 2;
		if ($mid > 600) $mid = 600;

		for ($i = 1;$i < count($lines) - 2;$i++)
		{
			$line = explode(" ", $lines[$i]);
			if (count($line) < 12) continue;
			$power = (float) $line[12];
			if ($power > $mid && $power > 650) break;
		}
		for (;$i > 1;$i--)
		{
			$line = explode(" ", $lines[$i]);
			$preline = explode(" ", $lines[$i - 1]);
			$power = (float) $line[12];
			$prepower = (float) $preline[12];
			if ($power < $prepower || $power < 500)
			{
				$runstart = $line[0];
				break;
			}
		}

		for ($i = count($lines) - 1;$i > 0;$i--)
		{
			$line = explode(" ", $lines[$i]);
			if (count($line) < 12) continue;
			$power = (float) $line[12];
			if ($power > $mid || $power > 500) break;
		}
		for (;$i < count($lines) - 1;$i++)
		{
			$line = explode(" ", $lines[$i]);
			$postline = explode(" ", $lines[$i + 1]);
			if (!isset($line[12]) || !isset($postpower[1])) continue;
			$power = (float) $line[12];
			$postpower = (float) $postline[12];
			if ($power < $postpower || $power < 550)
			{
				$runend = $line[0];
				break;
			}
		}

		if (isset($forcestarttime)) $runstart = $forcestarttime;
		if (isset($forceendtime)) $runend = $forceendtime;

		if (!isset($runstart) || !isset($runend))
		{
			echo "Cannot determine run duration in power log, skipping\n";
			continue;
		}

		if ($processingtime != -1)
		{
			echo "Time: Filter " . round($runend - $runstart, 2) . " Run $processingtime\n";
			if ($runend - $runstart > $processingtime && !isset($forceendtime)) $runend = $runstart + $processingtime;
		}

		$avgpower = $avgpower2 = $avgpower3 = $avgpower4 = $avgtime = $avgtime2 = $avgtime3 = $avgtime4 = 0;

		if ($starttime < $runstart - 0.15 * ($runend - $runstart)) $starttime = $runstart - 0.15 * ($runend - $runstart);
		if ($endtime > $runend + 0.15 * ($runend - $runstart)) $endtime = $runend + 0.15 * ($runend - $runstart);

		for ($i = 1;$i < count($lines);$i++)
		{
			$line = explode(" ", $lines[$i]);
			if (count($line) < 12) continue;
			$power = (float) $line[12];
			$time = (float) $line[0];
			$interval = (float) $line[3];
			if ($time >= $runstart && $time <= $runend)
			{
				$powerout[] = ($time - $runstart - $interval / 2) . " " . ($power * $interval) . " " . $interval;
				$avgpower += $power * $interval;
				$avgtime += $interval;
			}
			if ($time >= $runstart + 0.1 * ($runend - $runstart) && $time <= $runend)
			{
				$avgpower2 += $power * $interval;
				$avgtime2 += $interval;
			}
			if ($time >= $runstart + 0.7 * ($runend - $runstart) && $time <= $runend)
			{
				$avgpower3 += $power * $interval;
				$avgtime3 += $interval;
			}
			if ($time >= $runstart + 0.7 * ($runend - $runstart) && $time <= $runstart + 0.9 * ($runend - $runstart))
			{
				$avgpower4 += $power * $interval;
				$avgtime4 += $interval;
			}

			$plot[] = ($time - $starttime) . "\t$power";
		}
		if ($avgpower / $avgtime > 2000)
		{
			//$measurednodes = 2;
		}
		else
		{
			//$measurednodes = 1;
		}
		$avgpower /= $measurednodes;
		$avgpower2 /= $measurednodes;
		$avgpower3 /= $measurednodes;
		$avgpower4 /= $measurednodes;

		$plot = implode("\n", $plot);
		file_put_contents("plot$j.dat", $plot);

		$script .= '
		load "../pre.gnuplot"
		set style data lines
		set terminal pdf color enhanced font font_normal size 9.5in,5.5in dashed
		set output "' . str_replace(array(".mout", ".out"), "-power.pdf", $argv[$j]) . '"
		set xlabel "Time [s]"
		set ylabel "Power ' . ($measurednodes == 1 ? "per Node" : "per $measurednodes Nodes") . ' [W]"
		set key horizontal
		set xrange [0:' . ($endtime - $starttime) . ']
		set yrange [0:' . ($max * 1.17) . ']
		set arrow from ' . ($runstart - $starttime) . ',0 to ' . ($runstart - $starttime) . ',' . ($max) . 'nohead lw 4
		set arrow from ' . ($runend - $starttime) . ',0 to ' . ($runend - $starttime) . ',' . ($max) . 'nohead lw 4
		plot \'plot' . $j . '.dat\' using 1:2 notitle, \\
			-100 title "Result ' . $result . ' GFlop", \\
			' . $max . ' title "Max ' . round($max, 1) . ' W", \\
			-100 title "Avg ' . round($avgpower / $avgtime, 1) . ' W (' . round(1000 * $result / $avgpower * $avgtime / $numnodes, 1) . ' MFlop/W)"
			#-100 title "Green500 ' . round($avgpower2 / $avgtime2, 1) . ' W (' . round(1000 * $result / $avgpower2 * $avgtime2 / $numnodes, 1) . ' MFlop/W)"
		';
		echo "Efficiency $argv[$j]: Perf: $result GFlop, 0-100%: " . round($avgpower / $avgtime, 1) . ' W (' . round(1000 * $result / $avgpower * $avgtime / $numnodes, 1) . ' MFlop/W), 10-100%: ' . round($avgpower2 / $avgtime2, 1) . ' W (' . round(1000 * $result / $avgpower2 * $avgtime2 / $numnodes, 1) . ' MFlop/W), 70-100%: ' . round($avgpower3 / $avgtime3, 1) . ' W (' . round(1000 * $result / $avgpower3 * $avgtime3 / $numnodes, 1) . ' MFlop/W)' . ' MFlop/W), 70-90%: ' . round($avgpower4 / $avgtime4, 1) . ' W (' . round(1000 * $result / $avgpower4 * $avgtime4 / $numnodes, 1) . ' MFlop/W)' . "\n";
		
		$script .= '
		load "../pre.gnuplot"
		set style data lines
		set terminal pdf color enhanced font font_normal size 9.5in,5.5in dashed
		set output "' . str_replace(array(".mout", ".out"), "-power-70-90.pdf", $argv[$j]) . '"
		set xlabel "Time [s]"
		set ylabel "Power ' . ($measurednodes == 1 ? "per Node" : "per $measurednodes Nodes") . ' [W]"
		set key horizontal
		set xrange [0:' . ($endtime - $starttime) . ']
		set yrange [0:' . ($max * 1.17) . ']
		set arrow from ' . ($runstart - $starttime) . ',0 to ' . ($runstart - $starttime) . ',' . ($max) . 'nohead lw 4
		set arrow from ' . ($runend - $starttime) . ',0 to ' . ($runend - $starttime) . ',' . ($max) . 'nohead lw 4
		set arrow from ' . ($runstart - $starttime + 0.7 * ($runend - $runstart)) . ',0 to ' . ($runstart - $starttime + 0.7 * ($runend - $runstart)) . ',' . ($max) . 'nohead lw 4 lc rgb "green"
		set arrow from ' . ($runstart - $starttime + 0.9 * ($runend - $runstart)) . ',0 to ' . ($runstart - $starttime + 0.9 * ($runend - $runstart)) . ',' . ($max) . 'nohead lw 4 lc rgb "green"
		plot \'plot' . $j . '.dat\' using 1:2 notitle, \\
			-100 title "Result ' . $result . ' GFlop", \\
			' . $max . ' title "Max ' . round($max, 1) . ' W", \\
			-100 title "Avg (70%-90%) ' . round($avgpower4 / $avgtime4, 1) . ' W (' . round(1000 * $result / $avgpower4 * $avgtime4 / $numnodes, 1) . ' MFlop/W)"
			#-100 title "Green500 ' . round($avgpower2 / $avgtime2, 1) . ' W (' . round(1000 * $result / $avgpower2 * $avgtime2 / $numnodes, 1) . ' MFlop/W)"
		';

		
		$newname = str_replace(array(".mout", ".out", ".power", ".temp"), ".pfilter", $argv[$j]);
		$powerout[] = "EFF " . round(1000 * $result / $avgpower * $avgtime / $numnodes, 1);
		file_put_contents($newname, implode("\n", $powerout));
	}

	if (file_exists(str_replace(array(".mout", ".out", ".power", ".temp"), ".temp", $argv[$j])))
	{
		$file = file_get_contents(str_replace(array(".mout", ".out"), ".temp", $argv[$j]));
		$lines = explode("\n", str_replace("\r", '', $file));
		$time = -1;
		$adapter = -1;
		$plot = array();
		$entries = array();
		$peak = 0;
		$peaks = array();
		foreach ($lines as $line)
		{
			if ($line == '') continue;
			if (strpos($line, 'Adapter') === 0)
			{
				$tmp = explode(" ", $line, 3);
				$adapter = $tmp[1];
			}
			else if ($line{0} != ' ')
			{
				if ($time != -1)
				{
					$plot[$time] = $entries;
					$num_entries = count($entries);
					$entries = array();
				}
				while(strpos($line, "  ")) $line = str_replace("  ", " ", $line);
				$tmp = explode(" ", $line);
				if (is_numeric($tmp[1]{0}))
					$tmp = "$tmp[3] $tmp[1] $tmp[2] $tmp[5]";
				else
					$tmp = "$tmp[3] $tmp[2] $tmp[1] $tmp[5]";
				$time = strtotime($tmp);
				if ($time === false) echo "Error obtaining timestamp: $line $tmp\n";
			}
			else if (strpos($line, "Temperature") !== false)
			{
				while (strpos($line, "  ") !== false) $line = str_replace("  ", " ", $line);
				$tmp = explode(" ", $line);
				$entries[$adapter] = $tmp[5];
				if ($tmp[5] > $peak) $peak = $tmp[5];
				if (!isset($peaks[$adapter]) || $tmp[5] > $peaks[$adapter]) $peaks[$adapter] = $tmp[5];
			}
		}

		$lines = array();
		if (isset($starttime) && abs($time - $starttime) > 10000) unset($starttime);
		foreach ($plot as $time => $entries)
		{
			$lines[] = ($time - (isset($starttime) ? $starttime : (isset($starttimea) ? $starttimea : ($starttimea = $time)))) . "\t" . implode("\t", $entries) . "\t$time";
			if (isset($starttimea)) $endtimea = $time - $starttimea;
		}
		file_put_contents("plott$j.dat", implode("\n", $lines));

		$script .= '
		load "../pre.gnuplot"
		set style data lines
		set terminal pdf color enhanced font font_normal size 9.5in,5.5in solid
		set output "temperature_' . str_replace(array(".mout", ".out", ".power", ".temp"), ".pdf", $argv[$j]) . '"
		set xlabel "Time [s]"
		set ylabel "Temperature [°C]"
		set key horizontal
		set ytics 10
		set yrange[0:' . ($peak + 15) . ']
		' . (isset($starttime) ? ('set xrange [0:' . ($endtime - $starttime) . ']') : (isset($endtimea) && $endtimea ? "set xrange[0:$endtimea]" : '')) . '
		plot ';
		for ($k = 0;$k < $num_entries;$k++)
		{
			$script .= ' \'plott' . $j . '.dat\' using 1:' . (2 + $k) . ' title "GPU ' . "$k ($peaks[$k]" . '°C)",';
		}
		$script .= " $peak title \"Max Temperature $peak" . '°' . "C\" lc rgb \"#777777\"\n";

		if (isset($starttimea)) unset($starttimea);
	}	
}

file_put_contents("make.gnuplot", $script);

?>