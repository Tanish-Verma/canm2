# CANM2 Homework 2

This README gives a basic overview of what each folder/file is used for.

## Folder structure (from project root)

```

   hw2/
    ├── README.md
    ├── co2020_pdes_hw2.pdf
    ├── input1.in
    ├── input2.in
    ├── q1.c
    ├── q2.c
    ├── run_q1.sh
    ├── run_q2b.sh
    ├── run_q2c.sh
    ├── plotting_grids.py
    ├── plotting_numerical_analytical.py
    ├── plotting_q1.py
    ├── plotting_q2a.py
    ├── plotting_q2b.py
    ├── plotting_q2c.py
    ├── requirements.txt
    ├── output1/
    ├── output2/
    ├── plots/
```

## `hw2/` files

- `hw2/README.md`  
	This file.

- `hw2/co2020_pdes_hw2.pdf`  
	Assignment PDF.

- `hw2/requirements.txt`  
	Python dependencies for plotting scripts.

- `hw2/input1.in`  
	Input file used by `q1.c`.

- `hw2/input2.in`  
	Input file used by `q2.c`.

- `hw2/q1.c`  
	Main C solver for Q1 on stretched grid, writes solution and error data in `output1/`.

- `hw2/q2.c`  
	Main C solver for Q2 (stability/convergence runs), writes outputs in `output2/`.

- `hw2/run_q1.sh`  
	Helper script to compile and run Q1 cases.

- `hw2/run_q2b.sh`  
	Helper script for Q2(b) stability runs.

- `hw2/run_q2c.sh`  
	Helper script for Q2(c) convergence runs.

- `hw2/plotting_grids.py`  
	Plots uniform/stretched grids and spacing distributions.

- `hw2/plotting_numerical_analytical.py`  
	Simple comparison plot for numerical vs analytical solution on the finest grid.

- `hw2/plotting_q1.py`  
	Q1 convergence plotting from `output1/error_b*.dat` files.

- `hw2/plotting_q2a.py`  
	Von Neumann stability plotting for different `r` values.

- `hw2/plotting_q2b.py`  
	Q2(b) time-evolution visualization for stable/unstable cases.

- `hw2/plotting_q2c.py`  
	Q2(c) convergence and saturation behavior plot.

## `hw2/` subfolders

- `hw2/output1/`  
	Q1 outputs (`T_xy_*.dat`, `error*.dat`).

- `hw2/output2/`  
	Q2 outputs (time-step snapshots and error data).

- `hw2/plots/`  
	Generated plot images from all `plotting_*.py` scripts.

---
