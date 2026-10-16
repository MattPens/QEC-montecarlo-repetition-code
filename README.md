# Monte Carlo Simulation of a 3-Qubit Quantum Error-Correction Code

Numerical study of the 3-qubit repetition code under independent bit-flip
noise. The project is written in C++ as an incremental exercise in numerical
algorithms, Monte Carlo methods, stabilizer syndromes, decoding, and error
analysis.


## Goal

Estimate with Monte Carlo the logical error probability `P_L(p)` of an ideal
3-qubit repetition code and compare it with the exact analytical result.

Each physical qubit undergoes an independent bit flip with probability `p`.
The code measures the stabilizers `Z1Z2` and `Z2Z3`, identifies a syndrome,
applies a minimum-weight correction, and classifies the final state as either
a successful correction or a logical error.

## Model

For a physical error pattern `e = (e0,e1,e2)`, each `ei` is either `0`
(no bit flip) or `1` (bit flip). The error weight is

```text
w = e0 + e1 + e2.
```

The decoder corrects all patterns of weight 0 or 1. Patterns of weight 2 or 3
produce the residual pattern `111`, corresponding to a logical bit flip.

Therefore,

```math
P_L(p) = P(w=2) + P(w=3)
       = 3 p^2 (1-p) + p^3
       = 3 p^2 - 2 p^3.
```

For `p < 0.5`, the code improves reliability because `P_L(p) < p`.

## Monte Carlo procedure

For each trial:

1. Generate three independent Bernoulli errors with probability `p`.
2. Compute the two stabilizer syndrome bits.
3. Decode the syndrome and apply the correction.
4. Check whether the residual pattern is `111`.
5. Repeat `N` times and estimate `P_L = K_L / N`.

The estimated statistical uncertainty is

```math
sigma(P_L) = sqrt(P_L (1-P_L) / N).
```

## Validation

The program includes several checks:

- convergence of uniform random moments and Bernoulli frequencies;
- comparison of sampled error-weight probabilities with the binomial law;
- exhaustive verification of all eight 3-qubit error configurations;
- agreement between Monte Carlo and the analytical logical-error curve;
- convergence of the Monte Carlo uncertainty as `1/sqrt(N)`.

For example, at `p = 0.1` and `N = 100000`, a typical run gives

```text
P_L(Monte Carlo) = 0.02894 +/- 0.00053
P_L(theory)      = 0.02800
```

The estimate is statistically compatible with the exact result.

## Build and run

```bash
g++ -std=c++17 -Wall -Wextra -pedantic qec_montecarlo.cpp -o qec_progetto
./qec_progetto
```

The program generates `curva_PLogica.dat`, with columns

```text
p   P_L_MonteCarlo   uncertainty   P_L_theory
```

## Plot with Gnuplot

```gnuplot
set xlabel "Physical error probability p"
set ylabel "Logical error probability"
set xrange [0:1]
set yrange [0:1]
set grid
set key top left

plot "curva_PLogica.dat" using 1:2:3 with yerrorbars \
     title "Monte Carlo", \
     "curva_PLogica.dat" using 1:4 with lines lw 2 \
     title "Theory", \
     x with lines lw 2 title "Unencoded: P_L = p"
```

## Scope and limitations

This is an idealized QEC simulation: syndrome extraction and correction are
assumed perfect. The model includes only independent bit-flip noise, not
phase-flip, depolarizing noise, faulty gates, faulty measurements, ancilla
qubits, or a fault-tolerant circuit.

Possible extensions include phase-flip protection, Pauli/depolarizing
channels, larger stabilizer codes, and circuit-level fault-tolerant models.

## Author

Mattia Pennone 
MSc student in Complex Systems Physics - Quantum Computation, Communication and Machine
Learning, University of Turin.