#!/bin/bash

# -------------------------------------------------------
# Runs hc1d.c for 4 values of r using nx=50
# dx = 1/(nx-1) = 1/49 ~ 0.02041
# dt = r * dx^2
#   r=0.40 -> dt ~ 1.6639e-4
#   r=0.49 -> dt ~ 2.0383e-4
#   r=0.51 -> dt ~ 2.1227e-4
#   r=0.70 -> dt ~ 2.9121e-4
# -------------------------------------------------------

NX=50
XST=0.0
XEN=1.0
TST=0.0
TEN=0.05
OUTBASE="output2"

mkdir -p "$OUTBASE/r040" "$OUTBASE/r049" "$OUTBASE/r051" "$OUTBASE/r070"

cat > input2.in << INPUTEOF
$NX
$XST $XEN
$TST $TEN
INPUTEOF

run_case() {
    local dt=$1
    local folder=$2
    local outdir="$OUTBASE/$folder"

    echo "--- Running case: dt=$dt -> $outdir ---"
    sed -i "s/^\([[:space:]]*\)dt = [^;]*;/\1dt = $dt;/" q2.c
    gcc -o q2 q2.c -lm
    if [ $? -ne 0 ]; then echo "Compilation failed!"; exit 1; fi
    rm -f output2/output_*.dat
    rm -f "$outdir"/output_*.dat
    ./q2
    mv output2/output_*.dat "$outdir"/
    echo ""
}

run_case "1.6639e-4" "r040"
run_case "2.0383e-4" "r049"
run_case "2.1227e-4" "r051"
run_case "2.9121e-4" "r070"

echo "All runs done. Run: python plotting_q2b.py"