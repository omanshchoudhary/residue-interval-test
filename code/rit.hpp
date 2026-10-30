// Dependence tests for one linear equation:
//   a[0]*x[0] + a[1]*x[1] + ... = c,  with lo[k] <= x[k] <= hi[k]
// Every test returns true for "a dependence may exist" and false only when
// it has proved that no integer solution exists.
#pragma once

#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <vector>
using namespace std;

struct Problem {
    vector<long long> a, lo, hi;
    long long c = 0;
};

long long floor_div(long long a, long long b) {
    long long q = a / b;
    bool has_remainder = a % b != 0;
    bool negative = (a < 0) != (b < 0);
    if (has_remainder && negative) {
        q--;
    }
    return q;
}

long long ceil_div(long long a, long long b) {
    return -floor_div(-a, b);
}

// remainder that is never negative
long long mod(long long a, long long m) {
    return ((a % m) + m) % m;
}

// smallest and largest value of a*x when x is in [lo, hi]
long long term_min(long long a, long long lo, long long hi) {
    if (a > 0) {
        return a * lo;
    }
    return a * hi;
}

long long term_max(long long a, long long lo, long long hi) {
    if (a > 0) {
        return a * hi;
    }
    return a * lo;
}

// returns gcd(a, b) and sets x, y so that a*x + b*y = gcd(a, b)
long long extended_gcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    long long x1, y1;
    long long g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

// gcd of all coefficients must divide c
bool gcd_test(const Problem &p) {
    long long g = 0;
    for (long long value : p.a) {
        g = gcd(g, llabs(value));
    }
    if (g == 0) {
        return p.c == 0;
    }
    return p.c % g == 0;
}

// c must lie between the smallest and largest possible left side
bool banerjee_test(const Problem &p) {
    long long low = 0, high = 0;
    for (size_t i = 0; i < p.a.size(); i++) {
        low += term_min(p.a[i], p.lo[i], p.hi[i]);
        high += term_max(p.a[i], p.lo[i], p.hi[i]);
    }
    return low <= p.c && p.c <= high;
}

// the I test, written from secondary descriptions: divide by the gcd, and move
// a term into the interval once its coefficient is at most the interval length
bool i_test(const Problem &p) {
    vector<long long> a = p.a, lo = p.lo, hi = p.hi;
    long long low = p.c, high = p.c;

    while (true) {
        if (a.empty()) {
            return low <= 0 && 0 <= high;
        }

        long long g = 0;
        for (long long value : a) {
            g = gcd(g, llabs(value));
        }
        low = ceil_div(low, g);
        high = floor_div(high, g);
        if (low > high) {
            return false;
        }
        for (size_t i = 0; i < a.size(); i++) {
            a[i] /= g;
        }

        bool moved = false;
        for (size_t k = 0; k < a.size(); k++) {
            if (llabs(a[k]) <= high - low + 1) {
                low -= term_max(a[k], lo[k], hi[k]);
                high -= term_min(a[k], lo[k], hi[k]);
                a.erase(a.begin() + k);
                lo.erase(lo.begin() + k);
                hi.erase(hi.begin() + k);
                moved = true;
                break;
            }
        }

        // nothing can move, so fall back to a bounds check on what is left
        if (!moved) {
            long long sum_low = 0, sum_high = 0;
            for (size_t i = 0; i < a.size(); i++) {
                sum_low += term_min(a[i], lo[i], hi[i]);
                sum_high += term_max(a[i], lo[i], hi[i]);
            }
            return sum_low <= high && low <= sum_high;
        }
    }
}

// rule 1: a*x = c modulo g, so x keeps one remainder modulo g / gcd(a, g).
// moves new_lo up and new_hi down to the nearest allowed values.
// assumes coefficients are small enough that m*m fits in 64 bits.
void apply_residue_rule(long long a, long long c, long long g,
                        long long &new_lo, long long &new_hi) {
    if (g <= 1) {
        return;
    }
    long long inverse, unused;
    long long d = extended_gcd(mod(a, g), g, inverse, unused);
    long long m = g / d;
    if (m <= 1) {
        return;
    }
    long long reduced_c = mod(mod(c, g) / d, m);
    long long remainder = reduced_c * mod(inverse, m) % m;
    new_lo += mod(remainder - new_lo, m);
    new_hi -= mod(new_hi - remainder, m);
}

// rule 2: a*x must lie in [from, to], so divide and round inward
void apply_interval_rule(long long a, long long from, long long to,
                         long long &new_lo, long long &new_hi) {
    if (a > 0) {
        new_lo = max(new_lo, ceil_div(from, a));
        new_hi = min(new_hi, floor_div(to, a));
    } else {
        new_lo = max(new_lo, ceil_div(to, a));
        new_hi = min(new_hi, floor_div(from, a));
    }
}

// rule 3: terms not divisible by g must add up to c modulo g.
// false means no value in their range can do that.
bool block_rule_holds(const Problem &p, const vector<long long> &lo,
                      const vector<long long> &hi) {
    int n = p.a.size();

    // try each coefficient and each pairwise gcd as the divisor
    vector<long long> divisors;
    for (int i = 0; i < n; i++) {
        divisors.push_back(llabs(p.a[i]));
        for (int j = i + 1; j < n; j++) {
            divisors.push_back(gcd(llabs(p.a[i]), llabs(p.a[j])));
        }
    }
    sort(divisors.begin(), divisors.end());
    divisors.erase(unique(divisors.begin(), divisors.end()), divisors.end());

    for (long long g : divisors) {
        if (g < 2) {
            continue;
        }
        long long low = 0, high = 0;
        for (int i = 0; i < n; i++) {
            if (p.a[i] % g != 0) {
                low += term_min(p.a[i], lo[i], hi[i]);
                high += term_max(p.a[i], lo[i], hi[i]);
            }
        }
        long long first_match = low + mod(p.c - low, g);
        if (first_match > high) {
            return false;
        }
    }
    return true;
}

// shrink the ranges with the three rules until nothing changes.
// an empty range means no solution exists.
bool residue_interval_test(const Problem &p, int *rounds_used = nullptr,
                           int round_limit = 1000) {
    if (!gcd_test(p)) {
        return false;
    }
    int n = p.a.size();
    vector<long long> lo = p.lo, hi = p.hi;

    // prefix and suffix gcds give "gcd of all but k" without recomputing it
    vector<long long> prefix(n + 1, 0), suffix(n + 2, 0);
    for (int i = 0; i < n; i++) {
        prefix[i + 1] = gcd(prefix[i], llabs(p.a[i]));
    }
    for (int i = n - 1; i >= 0; i--) {
        suffix[i] = gcd(suffix[i + 1], llabs(p.a[i]));
    }

    for (int round = 1; round <= round_limit; round++) {
        if (rounds_used) {
            *rounds_used = round;
        }
        bool changed = false;

        long long low = 0, high = 0;
        for (int i = 0; i < n; i++) {
            low += term_min(p.a[i], lo[i], hi[i]);
            high += term_max(p.a[i], lo[i], hi[i]);
        }

        for (int k = 0; k < n; k++) {
            long long a = p.a[k];

            // take term k out so low and high cover the other terms only
            low -= term_min(a, lo[k], hi[k]);
            high -= term_max(a, lo[k], hi[k]);

            long long new_lo = lo[k], new_hi = hi[k];
            long long g = gcd(prefix[k], suffix[k + 1]);
            apply_residue_rule(a, p.c, g, new_lo, new_hi);
            apply_interval_rule(a, p.c - high, p.c - low, new_lo, new_hi);

            if (new_lo > new_hi) {
                return false;
            }
            if (new_lo != lo[k] || new_hi != hi[k]) {
                changed = true;
            }
            lo[k] = new_lo;
            hi[k] = new_hi;

            low += term_min(a, lo[k], hi[k]);
            high += term_max(a, lo[k], hi[k]);
        }

        if (!block_rule_holds(p, lo, hi)) {
            return false;
        }
        if (!changed) {
            break;
        }
    }
    return true;
}

// the exact answer: track every sum reachable so far
bool has_solution(const Problem &p) {
    int n = p.a.size();
    long long low = 0, high = 0;
    for (int i = 0; i < n; i++) {
        low += term_min(p.a[i], p.lo[i], p.hi[i]);
        high += term_max(p.a[i], p.lo[i], p.hi[i]);
    }
    if (p.c < low || p.c > high) {
        return false;
    }

    // index s stands for the sum s + low
    long long width = high - low + 1;
    vector<char> reachable(width, 0), next(width);
    reachable[-low] = 1;
    for (int i = 0; i < n; i++) {
        fill(next.begin(), next.end(), 0);
        for (long long s = 0; s < width; s++) {
            if (!reachable[s]) {
                continue;
            }
            for (long long x = p.lo[i]; x <= p.hi[i]; x++) {
                long long t = s + p.a[i] * x;
                if (t >= 0 && t < width) {
                    next[t] = 1;
                }
            }
        }
        swap(reachable, next);
    }
    return reachable[p.c - low];
}
