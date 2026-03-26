#!/bin/bash

# Domain boundaries
XST=0.0
XEN=1.0
YST=0.0
YEN=0.5

# Grid sizes to test
grid_sizes=(5 9 17 33 65 129 257 513)

# Set beta values here once; script and output naming follow automatically
BETAS=(1.5 3.0)

# Create output directory if it doesn't exist
mkdir -p output1

# -------------------------------------------------------
# Function: run all grid sizes for a given compiled binary
# and save the resulting error.dat to a named file
# -------------------------------------------------------
run_all_grids() {
    local output_file=$1

    # Clean any previous error.dat
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

generated_files=()

for beta in "${BETAS[@]}"; do
    beta_tag=$(awk -v b="$beta" 'BEGIN { printf "%03d", int((b * 100) + 0.5) }')

    echo "=========================================="
    echo "Compiling with BETA = $beta ..."
    echo "=========================================="
    sed -i "s/#define BETA .*/#define BETA $beta/" q1.c
    gcc -o a.out q1.c -lm
    if [ $? -ne 0 ]; then echo "Compilation failed!"; exit 1; fi
    echo "Compilation successful!"
    echo ""

    output_file="error_b${beta_tag}.dat"
    run_all_grids "$output_file"
    generated_files+=("output1/$output_file -- errors for beta=$beta")
done

# -------------------------------------------------------
echo "=========================================="
echo "All runs complete!"
for line in "${generated_files[@]}"; do
    echo "  $line"
done
echo "Run: python plotting_q1_convergence.py"
echo "=========================================="