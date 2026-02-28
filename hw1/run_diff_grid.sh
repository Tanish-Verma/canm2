#!/bin/bash

# Script to run q4 for different grid sizes
# Input format: nx ny
#               xst xen
#               yst yen

# Compile the code first
echo "Compiling q4.c..."
gcc -o q4 q4.c -lm
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Compilation successful!"
echo ""
# Clean the error.dat file in output folder
if [ -f output/error.dat ]; then
    rm output/error.dat
fi
echo ""
# Define grid sizes to test
grid_sizes=(5 9 17 33 65 129 257 513)

# Domain boundaries
XST=0.0
XEN=1.0
YST=0.0
YEN=0.5

# Run for each grid size
for nx in "${grid_sizes[@]}"; do
    ny=$((nx/2 + 1))
    
    echo "=========================================="
    echo "Running with grid size: $nx x $ny"
    echo "=========================================="
    echo ""
    
    # Create input.in
    cat > input.in << EOF
$nx $ny
$XST $XEN
$YST $YEN
EOF
    
    echo "Created input.in with nx=$nx, ny=$ny"
    
    # Run the program
    ./q4
    
    echo ""
done

echo "=========================================="
echo "All grid sizes completed!"
echo "=========================================="


