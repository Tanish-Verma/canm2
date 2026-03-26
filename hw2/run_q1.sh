#!/bin/bash

XST=0.0
XEN=1.0
YST=0.0
YEN=0.5

grid_sizes=(5 9 17 33 65 129 257 513)
BETAS=(1.5 3.0)

mkdir -p output1

run_all_grids() {
    rm -f output1/error.dat
    for nx in "${grid_sizes[@]}"; do
        ny=$((nx / 2 + 1))

        echo "=========================================="
        echo "Running with grid size: $nx x $ny"
        echo "=========================================="

        cat > input1.in << EOF
$nx $ny
$XST $XEN
$YST $YEN
EOF
        ./a.out

        echo ""
    done

    # Save collected errors for this beta
    cp output1/error.dat "output1/$output_file"
    echo ">> Saved errors to output1/$output_file"
    echo ""
}

for beta in "${BETAS[@]}"; do
    beta_tag=$(awk -v b="$beta" 'BEGIN { printf "%03d", int((b * 100) + 0.5) }')
    echo "Compiling with BETA = $beta ..."
    sed -i "s/#define BETA .*/#define BETA $beta/" q1.c
    gcc -o a.out q1.c -lm
    if [ $? -ne 0 ]; then echo "Compilation failed!"; exit 1; fi
    run_all_grids "error_b${beta_tag}.dat"
done

echo "Done! Run: python plotting_q1_convergence.py"