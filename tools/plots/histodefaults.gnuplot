set output GPVAL_OUTPUT
ntics = GPVAL_DATA_X_MAX - 0.5
nboxes_tmp = nboxes == 1 ? 1 : nboxes - 1
set xtics nomirror scale 0
set xrange [-1.+1./(2.*nboxes_tmp):ntics - 1./(2.*nboxes_tmp)]
