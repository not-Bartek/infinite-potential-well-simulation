# Quantum Particle Simulation

A C++17 simulation of a particle in a one-dimensional infinite potential well. The
program combines user-defined stationary states and calculates the probability
density over time.

## Requirements

- GNU Make
- A C++17 compiler (`g++`)

macOS, Linux, and Windows are supported. On Windows, the program can be
compiled directly with a C++17 compiler such as `g++`. The Makefile requires a
POSIX-compatible shell, such as WSL, MSYS2, or Git Bash.

## Build and run

```bash
make
make run
```

The program first asks for the number of states. For each state, enter:

- `n` - stationary-state number
- `phi` - phase
- `c` - amplitude

The amplitudes are normalized before the simulation starts.

## Simulation parameters

Default values are:

```text
well length (L): 10.0
spatial intervals: 500
time step (dt): 0.01
maximum time (tMax): 10.0
particle mass: 9.110e-31
```

Override the first four values from the command line:

```bash
./program <L> <intervals> <dt> <tMax>
```

Example:

```bash
./program 0.1 500 0.01 100
```

### Optional: compile without Make

```bash
g++ -std=c++17 -Wall -Wextra -pedantic main.cpp -o program
./program 0.1 500 0.01 100
```

## Output

The simulation writes probability-density data to `simulation.csv`. Each row
contains the time followed by the density at every spatial interval. Existing
CSV files in the repository are example simulation results.

## Visualization with `visualize.py`

`visualize.py` reads the simulation CSV file and displays an animated plot of
the probability density over the well. It requires Python with `pandas`,
`numpy`, and `matplotlib` installed.

```bash
python3 visualize.py
python3 visualize.py simulation.csv 0.1
```

The optional arguments are the CSV path and the well length `L`. Use the same
`L` value as in the simulation command.

## Cleaning

```bash
make clean
```

## License

This project is licensed under the MIT License.

Copyright (c) 2026

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
