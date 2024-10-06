# Numerical Analysis Homework #2

This project implements the bisection method, Newton's method, and the secant method in C++ to solve nonlinear equations.

## Files

- `Function.hpp`: Declaration of the `Function` class.
- `EquationSolver.hpp`: Declaration and implementation of the equation solver classes.
- `main.cpp`: Main program with test cases.
- `report.tex`: LaTeX source code for the report.
- `Makefile`: Build and run instructions.
- `README.md`: Project description.
- `.gitignore`: Files to ignore in the repository.

## Build and Run

To compile the program:

```bash
make
```

To run the program and save the output to output.txt:

```bash
make run
```

To compile the report into a PDF:

```bash
make report
```

To clean up generated files:

```bash
make clean
```

## Requirements
C++ compiler (e.g., g++)
LaTeX distribution (e.g., TeX Live, MiKTeX) for compiling the report.
## Description
The program solves various nonlinear equations using the implemented numerical methods and outputs the results. The report provides a detailed explanation and analysis of the results.