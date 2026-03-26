#!/bin/bash

XST=0.0
XEN=1.0
TST=0.001
TEN=0.002

grid_sizes=(10 20 40 80 160 320 640 1000)

mkdir -p output2
rm -f output2/error2.dat

echo "Compiling q2.c..."
sed -i "s/^\([[:space:]]*\)dt = [^;]*;/\1dt = 1.25e-7;/" q2.c
gcc -o q2 q2.c -lm
if [ $? -ne 0 ]; then echo "Compilation failed!"; exit 1; fi

for nx in "${grid_sizes[@]}"; do
    cat > input2.in << INPUTEOF
$nx
$XST $XEN
$TST $TEN
INPUTEOF
    ./q2
done

echo "Done! Run: python plotting_q2c.py"