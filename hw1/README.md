# CANM2 Homework 1

This README gives a basic overview of what each folder/file is used for.

## Folder structure (from project root)

```

   hw1/
    ├── README.md
    ├── co2020_pdes_hw1.pdf
    ├── input.in
    ├── q3.py
    ├── q4.c
    ├── q4_plotting.py
    ├── run_diff_grid.sh
    ├── output/
    ├── plots/
```

## `hw1/` files

- `requirements.txt`  
	Python dependencies (mainly numpy/matplotlib + notebook tooling).


- `hw1/README.md`  
	This file.

- `hw1/co2020_pdes_hw1.pdf`  
	Assignment PDF .

- `hw1/input.in`  
	Input file read by `q4.c`.

- `hw1/q3.py`  
	Python script for analytical solution plotting (contour + x/y profiles).

- `hw1/q4.c`  
	Main C solver made using the skeleton code, that computes numerical solution, exact solution, and L2 error.

- `hw1/q4_plotting.py`  
	Python post-processing script to plot:
	- numerical contour from `output/T_xy_*.dat`
	- error convergence from `output/error.dat`
	- profile comparisons across grid sizes.

- `hw1/run_diff_grid.sh`  
	Helper script to compile `q4.c`, run it for multiple grid sizes, and regenerate `output/error.dat`.

## `hw1/` subfolders

- `hw1/output/`  
	Numerical output data files (`T_xy_*.dat`) and `error.dat`.

- `hw1/plots/`  
	Generated plot images from `q3.py` and `q4_plotting.py`.

---
