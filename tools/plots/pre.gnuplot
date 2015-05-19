reset

font_normal = "CMU Sans Serif,10"
font_small = "CMU Sans Serif,8"
font_smaller = "CMU Sans Serif,7"
font_smallest = "CMU Sans Serif,6"

set terminal pdf color enhanced font font_normal size 5in,3in solid
set style histogram clustered gap 1
set style fill solid border -1
set grid y
set bars front
set tics front
set grid back
set bars 2
set mytics 5
set key right top

set autoscale y
set autoscale x
set ytics auto
set xtic auto

user_color1 = "#B03030"
#176 48 48
user_color2 = "#00A000"
#0 160 0
user_color3 = "#0000C0"
#0 0 192
user_color4 = "#9400D3"
#148 0 211
user_color5 = "#19BBBF"
#25 187 191
user_color6 = "#F25900"
user_color7 = "black"
user_color8 = "#7F7F7F"
user_color9 = "#FFD700"
user_color10 = "#07F707"
user_color11 = "#07F7F7"
user_color12 = "#F08080"

#user_color1 = "#000000"
#user_color2 = "#FFFFFF"
#user_color3 = "#555555"
#user_color4 = "#AAAAAA"

user_bgcolor1 = "#CA7575"
user_bgcolor2 = "#55BF55"
user_bgcolor3 = "#5555D5"
user_bgcolor4 = "#B755E1"
user_bgcolor5 = "#65D1D4"
user_bgcolor6 = "#F69055"
user_bgcolor7 = "#555555"
user_bgcolor8 = "#A9A9A9"
user_bgcolor9 = "#FFE455"
user_bgcolor10 = "#59F959"
user_bgcolor11 = "#59F9F9"
user_bgcolor12 = "#F5AAAAs"


my_linewidth=0.5
my_pointsize=0.75

set pointsize my_pointsize
nboxes=1
valpos(a)=($0-1+(a-(nboxes-1.)/2.)/(nboxes+1.))

set style line 1 lw my_linewidth lt rgb user_color1 pt 1
set style line 2 lw my_linewidth lt rgb user_color2 pt 2
set style line 3 lw my_linewidth lt rgb user_color3 pt 3
set style line 4 lw my_linewidth lt rgb user_color4 pt 4
set style line 5 lw my_linewidth lt rgb user_color5 pt 6
set style line 6 lw my_linewidth lt rgb user_color6 pt 7
set style line 7 lw my_linewidth lt rgb user_color7 pt 8
set style line 8 lw my_linewidth lt rgb user_color8 pt 9
set style line 9 lw my_linewidth lt rgb user_color9 pt 10
set style line 10 lw my_linewidth lt rgb user_color10 pt 11
set style line 11 lw my_linewidth lt rgb user_color11 pt 5
set style line 12 lw my_linewidth lt rgb user_color12 pt 12
set style line 13 lw my_linewidth lt rgb user_color1 pt 13
set style line 14 lw my_linewidth lt rgb user_color2 pt 14
set style line 15 lw my_linewidth lt rgb user_color3 pt 15
set style line 16 lw my_linewidth lt rgb user_color4 pt 16
set style line 17 lw my_linewidth lt rgb user_color5 pt 17
set style line 18 lw my_linewidth lt rgb user_color6 pt 18
set style line 19 lw my_linewidth lt rgb user_color7 pt 19
set style line 20 lw my_linewidth lt rgb user_color8 pt 20
set style line 21 lw my_linewidth lt rgb user_color9 pt 21
set style line 22 lw my_linewidth lt rgb user_color10 pt 22
set style line 23 lw my_linewidth lt rgb user_color11 pt 23
set style line 24 lw my_linewidth lt rgb user_color12 pt 24
set style increment user

unset xlabel
unset colorbox
set ylabel "MISSING"

set grid lt 0 lw 1 lc rgb '#8f8f8f', lt 0 lw 1 lc rgb '#8f8f8f'

set tmargin 0.5
