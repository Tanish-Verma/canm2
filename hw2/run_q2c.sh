#!/bin/bash

# -------------------------------------------------------
# Convergence study for Q2c
# Fixed: tst=0.001, ten=0.002, dt=1.25e-7 (already in hc1d.c)
# Runs for nx = 10, 20, 40, 80, 160, 320, 640, 1000
# -------------------------------------------------------

XST=0.0
XEN=1.0
TST=0.001
TEN=0.002

grid_sizes=(10 20 40 80 160 320 640 1000)

mkdir -p output2

# Clear previous error file
rm -f output2/error2.dat

echo "Compiling q2.c..."
sed -i "s/^\([[:space:]]*\)dt = [^;]*;/\1dt = 1.25e-7;/" q2.c
gcc -o q2 q2.c -lm
if [ $? -ne 0 ]; then echo "Compilation failed!"; exit 1; fi
echo "Compilation successful!"
echo ""

for nx in "${grid_sizes[@]}"; do
    echo "=========================================="
    echo "Running nx = $nx"
    echo "=========================================="

    cat > input2.in << INPUTEOF
$nx
$XST $XEN
$TST $TEN
INPUTEOF

    ./q2
    echo ""
done

echo "=========================================="
echo "All runs done!"
echo "Errors saved to output2/error2.dat"
echo "Run: python plot_q2c.py"
echo "=========================================="