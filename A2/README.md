# CANM2 Assignment 2

This README gives a basic overview of what each folder/file is used for.

## Folder structure (from project root)

``` text
canm2/
└── A2/
    ├── README.md
    ├── Q1/
    ├── Q2/
    ├── Q3/
    ├── Q4/
    └── Q5/
```

## `A2/` files

- `A2/README.md`  
This file.

## `A2/` subfolders

- `A2/Q1/`  
Q1 implementation and outputs.
- `CO24BTECH11023_A02Q01.cpp`: C++ implementation.
- `CO24BTECH11023_A02Q01.py`: Python implementation.
- `CO24BTECH11023_A02Q01`: Compiled executable (generated from C++ source).
- `data/`: Generated numerical output files.
- `plots/`: Generated plots.

- `A2/Q2/`  
Q2 implementation and outputs.
- `CO24BTECH11023_A02Q02.cpp`: C++ implementation.
- `CO24BTECH11023_A02Q02.py`: Python implementation.
- `CO24BTECH11023_A02Q02`: Compiled executable (generated from C++ source).
- `data/`: Generated numerical output files.
- `plots/`: Generated plots.

- `A2/Q3/`  
Q3 implementation and outputs.
- `CO24BTECH11023_A02Q03.cpp`: C++ implementation.
- `CO24BTECH11023_A02Q03.py`: Python implementation.
- `CO24BTECH11023_A02Q03`: Compiled executable (generated from C++ source).
- `data/`: Generated numerical output files.
- `plots/`: Generated plots.

- `A2/Q4/`  
Q4 implementation and outputs.
- `CO24BTECH11023_A02Q04.cpp`: C++ implementation.
- `CO24BTECH11023_A02Q04.py`: Python implementation.
- `CO24BTECH11023_A02Q04`: Compiled executable (generated from C++ source).
- `data/`: Generated numerical output files.
- `plots/`: Generated plots.

- `A2/Q5/`  
Q5 implementation and outputs.
- `CO24BTECH11023_A02Q05.cpp`: C++ implementation.
- `CO24BTECH1103_A02Q05.py`: Python implementation.
- `CO24BTECH11023_A02Q05`: Compiled executable (generated from C++ source).
- `data/`: Generated numerical output files.
- `plots/`: Generated plots.

## How to run all solutions (quick check)

From project root (`canm2/`), use the commands below to compile and run all C++ and Python files:

```bash
for q in 1 2 3 4 5; do
echo "Running Q${q} C++"
(
    cd "A2/Q${q}" || exit 1
    src="CO24BTECH11023_A02Q0${q}.cpp"
    exe="CO24BTECH11023_A02Q0${q}"
    g++ -std=c++20 -O2 "$src" -o "$exe" && ./$exe
)

done

for q in 1 2 3 4 5; do
echo "Running Q${q} Python"
(
    cd "Q${q}" || exit 1
    python3 "CO24BTECH11023_A02Q0${q}.py"
)
done
```

## Notes

- Run each program from A2 folder, so that folder paths work
- Python dependencies (for plotting or post-processing) may require `numpy` and `matplotlib` depending on the script.
