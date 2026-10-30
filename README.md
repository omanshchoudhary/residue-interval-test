# Residue-Interval Test

A small data dependence test for loop parallelization. It decides whether an equation

    a1*x1 + a2*x2 + ... + an*xn = c,   Lk <= xk <= Uk,   every xk an integer

can have a solution, using only gcd calculations and interval arithmetic. If it proves that no solution exists, two array accesses in a loop can never touch the same element.

The GCD test ignores loop bounds, and the Banerjee test ignores that loop counters are integers. This test applies both facts together with three rules:

1. **Residue rule:** limit a variable to one remainder, using the gcd of the other coefficients.
2. **Interval rule:** narrow a variable's range using the range of the other terms.
3. **Block rule:** for a divisor g, check that the terms g does not divide can add up to c modulo g.

This is course work, not new research. Related ideas appear in the I test, the Omega test, the Power test and Granger's congruence analysis.

## Files

| File | Purpose |
|---|---|
| `code/rit.hpp` | The test, the GCD, Banerjee and I tests, and an exact checker |
| `code/rit_eval.cpp` | The evaluation program |
| `code/loop_outer.c` | Speed comparison on a flattened-array loop |
| `code/results.txt` | Saved output of the evaluation |
| `code/bench.txt` | Saved timings of the speed comparison |

## Build

Requires g++ and gcc with C++20 and OpenMP support.

    cd code
    g++ -O2 -std=c++20 rit_eval.cpp -o rit_eval

## Reproduce the results

Arguments are the number of variables, largest coefficient and largest upper bound. Random runs also take the number of trials and use a fixed seed.

Exhaustive table:

    ./rit_eval exhaustive 1 8 8
    ./rit_eval exhaustive 2 8 6
    ./rit_eval exhaustive 3 6 4
    ./rit_eval exhaustive 4 3 3

Random problems:

    ./rit_eval random 3 20 30 40000
    ./rit_eval random 3 100 100 8000
    ./rit_eval random 4 20 20 40000
    ./rit_eval random 4 50 40 8000
    ./rit_eval random 5 12 10 40000
    ./rit_eval random 6 8 8 40000

Worked examples:

    ./rit_eval single 3 0 10 5 0 10 7
    ./rit_eval single 1000 0 1999 1 0 499 -1000 0 1999 -1 0 499 500

Each variable takes three values (coefficient, lower bound, upper bound), followed by c.

Timing, which varies by machine:

    ./rit_eval timing

Speed comparison, also machine-dependent:

    gcc -O3 -fopenmp -DSERIAL loop_outer.c -o serial
    gcc -O3 -fopenmp -DSERIAL -ftree-parallelize-loops=4 loop_outer.c -o autopar
    gcc -O3 -fopenmp loop_outer.c -o outer
    ./serial && ./autopar && ./outer

GCC's report on the outer loop:

    gcc -O3 -fopenmp -DSERIAL -ftree-parallelize-loops=4 -c loop_outer.c -fdump-tree-parloops-details -o /dev/null
    grep FAILED loop_outer.c.*parloops*

The saved results came from GCC 15.2 on an Intel Core i7-10610U. A different GCC version may parallelize this loop differently.

## Limits

The test assumes linear subscripts, known constant bounds and no aliasing. It checks one equation at a time and gives no direction or distance information. It can answer "maybe" when no solution exists. The evaluation problems are synthetic, so the results do not predict behaviour on typical programs.
