// Evaluation program for the case study.
// Build: g++ -O2 -std=c++20 rit_eval.cpp -o rit_eval
#include "rit.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>

// counts how each test did over a set of problems
struct Tally {
    long long total = 0, dependent = 0, independent = 0;
    long long gcd_proved = 0, banerjee_proved = 0, both_proved = 0;
    long long i_proved = 0, rit_proved = 0, i_not_rit = 0, wrong = 0;
    int max_rounds = 0;
};

void add_problem(Tally &t, const Problem &p) {
    bool truth = has_solution(p);
    bool gcd_maybe = gcd_test(p);
    bool banerjee_maybe = banerjee_test(p);
    bool i_maybe = i_test(p);
    int rounds = 0;
    bool rit_maybe = residue_interval_test(p, &rounds);

    t.max_rounds = max(t.max_rounds, rounds);
    t.total++;

    if (truth) {
        t.dependent++;
        // a test that says independent here would be a real bug
        if (!gcd_maybe || !banerjee_maybe || !i_maybe || !rit_maybe) {
            t.wrong++;
        }
        return;
    }

    t.independent++;
    if (!gcd_maybe) t.gcd_proved++;
    if (!banerjee_maybe) t.banerjee_proved++;
    if (!gcd_maybe || !banerjee_maybe) t.both_proved++;
    if (!i_maybe) t.i_proved++;
    if (!rit_maybe) t.rit_proved++;
    if (!i_maybe && rit_maybe) t.i_not_rit++;
}

void print_tally(const Tally &t) {
    double percent = t.independent ? 100.0 / t.independent : 0;
    printf("total=%lld dependent=%lld independent=%lld\n", t.total, t.dependent, t.independent);
    printf("proved independent: GCD=%lld Banerjee=%lld GCD+Banerjee=%lld (%.2f%%) "
           "I=%lld (%.2f%%) RIT=%lld (%.2f%%)\n",
           t.gcd_proved, t.banerjee_proved, t.both_proved, t.both_proved * percent,
           t.i_proved, t.i_proved * percent, t.rit_proved, t.rit_proved * percent);
    printf("I proved but RIT missed=%lld  wrong answers=%lld  max rounds=%d\n",
           t.i_not_rit, t.wrong, t.max_rounds);
}

// every coefficient choice, every upper bound, every constant in range + 2 each side
void run_exhaustive(int n, long long max_coef, long long max_upper) {
    Tally t;
    bool miss_shown = false;

    vector<long long> values;
    for (long long v = -max_coef; v <= max_coef; v++) {
        if (v != 0) {
            values.push_back(v);
        }
    }
    long long value_count = values.size();

    long long coef_cases = 1, upper_cases = 1;
    for (int i = 0; i < n; i++) {
        coef_cases *= value_count;
        upper_cases *= max_upper;
    }

    Problem p;
    p.a.resize(n);
    p.lo.assign(n, 0);
    p.hi.resize(n);

    for (long long coef_case = 0; coef_case < coef_cases; coef_case++) {
        // read coef_case as a number in base value_count
        long long digits = coef_case;
        for (int i = 0; i < n; i++) {
            p.a[i] = values[digits % value_count];
            digits /= value_count;
        }

        for (long long upper_case = 0; upper_case < upper_cases; upper_case++) {
            digits = upper_case;
            for (int i = 0; i < n; i++) {
                p.hi[i] = 1 + digits % max_upper;
                digits /= max_upper;
            }

            long long low = 0, high = 0;
            for (int i = 0; i < n; i++) {
                low += term_min(p.a[i], 0, p.hi[i]);
                high += term_max(p.a[i], 0, p.hi[i]);
            }

            for (p.c = low - 2; p.c <= high + 2; p.c++) {
                add_problem(t, p);
                if (!miss_shown && !has_solution(p) && residue_interval_test(p)) {
                    miss_shown = true;
                    printf("first miss:");
                    for (int i = 0; i < n; i++) {
                        printf(" %+lldx%d[0,%lld]", p.a[i], i + 1, p.hi[i]);
                    }
                    printf(" = %lld\n", p.c);
                }
            }
        }
    }
    print_tally(t);
}

Problem random_problem(mt19937_64 &rng, int n, long long max_coef, long long max_upper) {
    Problem p;
    long long low = 0, high = 0;
    for (int i = 0; i < n; i++) {
        long long value = 1 + rng() % max_coef;
        if (rng() & 1) {
            value = -value;
        }
        long long upper = 1 + rng() % max_upper;
        p.a.push_back(value);
        p.lo.push_back(0);
        p.hi.push_back(upper);
        low += term_min(value, 0, upper);
        high += term_max(value, 0, upper);
    }

    // half the constants fall near the edge, where gaps in the sums live
    long long width = high - low + 1;
    if (rng() & 1) {
        width = min(width, 2 * max_coef + 1);
    }
    p.c = low + (long long)(rng() % width);
    return p;
}

void run_random(int n, long long max_coef, long long max_upper, long long trials) {
    mt19937_64 rng(12345 + n * 1000 + max_coef);
    Tally t;
    for (long long i = 0; i < trials; i++) {
        add_problem(t, random_problem(rng, n, max_coef, max_upper));
    }
    print_tally(t);
}

void run_single(int argc, char **argv) {
    Problem p;
    int n = (argc - 3) / 3;
    for (int i = 0; i < n; i++) {
        p.a.push_back(atoll(argv[2 + 3 * i]));
        p.lo.push_back(atoll(argv[3 + 3 * i]));
        p.hi.push_back(atoll(argv[4 + 3 * i]));
    }
    p.c = atoll(argv[argc - 1]);

    printf("GCD=%s Banerjee=%s I=%s RIT=%s truth=%s\n",
           gcd_test(p) ? "maybe" : "independent",
           banerjee_test(p) ? "maybe" : "independent",
           i_test(p) ? "maybe" : "independent",
           residue_interval_test(p) ? "maybe" : "independent",
           has_solution(p) ? "dependent" : "independent");
}

// wrappers so every test fits the same function pointer type
bool gcd_and_banerjee(const Problem &p) {
    return gcd_test(p) && banerjee_test(p);
}

bool rit_with_defaults(const Problem &p) {
    return residue_interval_test(p);
}

void time_test(const char *name, bool (*test)(const Problem &),
               const vector<Problem> &problems, size_t count) {
    auto start = chrono::steady_clock::now();
    long maybe_count = 0;
    for (size_t i = 0; i < count; i++) {
        maybe_count += test(problems[i]);
    }
    auto end = chrono::steady_clock::now();
    double ns = chrono::duration<double, nano>(end - start).count();
    printf("%-14s %10.0f ns per query (%ld)\n", name, ns / count, maybe_count);
}

void run_timing() {
    mt19937_64 rng(7);
    vector<Problem> problems;
    for (int i = 0; i < 200000; i++) {
        int n = 2 + rng() % 5;
        problems.push_back(random_problem(rng, n, 50, 100));
    }

    // three passes, so the last one is not skewed by a cold start
    for (int pass = 0; pass < 3; pass++) {
        time_test("GCD+Banerjee", gcd_and_banerjee, problems, problems.size());
        time_test("I test", i_test, problems, problems.size());
        time_test("RIT", rit_with_defaults, problems, problems.size());
        time_test("exact sums", has_solution, problems, 2000);
    }
}

int usage() {
    puts("usage: rit_eval exhaustive <vars> <max_coef> <max_upper>");
    puts("       rit_eval random <vars> <max_coef> <max_upper> <trials>");
    puts("       rit_eval single a1 lo1 hi1 a2 lo2 hi2 ... c");
    puts("       rit_eval timing");
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        return usage();
    }
    if (!strcmp(argv[1], "exhaustive") && argc == 5) {
        run_exhaustive(atoi(argv[2]), atoll(argv[3]), atoll(argv[4]));
    } else if (!strcmp(argv[1], "random") && argc == 6) {
        run_random(atoi(argv[2]), atoll(argv[3]), atoll(argv[4]), atoll(argv[5]));
    } else if (!strcmp(argv[1], "single") && argc >= 6 && (argc - 3) % 3 == 0) {
        run_single(argc, argv);
    } else if (!strcmp(argv[1], "timing")) {
        run_timing();
    } else {
        return usage();
    }
    return 0;
}
