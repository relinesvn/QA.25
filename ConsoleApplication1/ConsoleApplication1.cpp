// ConsoleApplication1.cpp
// Моделювання флуктуацій температури перегрітої пари (пароперегрівач)
// Версія для UNIT-тестів: основні функції винесені у model.h/model.cpp.
// Цей файл містить лише CLI/файловий ввід-вивід та "склеювання" моделі.

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <string>
#include <chrono>

#include "model.h"   // N, svertka(), normal12(), buildImpulse(), computeMeanSigma()

// ------------------------- CLI parsing -------------------------
// Usage:
//   ConsoleApplication1.exe <NN> [--seed=<uint>] [--no-file]
//
// Якщо <NN> не задано — запитати з консолі.
static bool parse_inputs(int argc, char* argv[], long& NN, bool& noFile, bool& hasSeed, uint32_t& seed)
{
    NN = 0;
    noFile = false;
    hasSeed = false;
    seed = 0;

    // 1) First positional arg: NN (if present and doesn't start with '-')
    if (argc >= 2 && argv[1] && argv[1][0] != '-') {
        char* endp = nullptr;
        long v = std::strtol(argv[1], &endp, 10);
        if (!endp || *endp != '\0') return false;
        NN = v;
    }

    // 2) Options
    for (int i = 1; i < argc; ++i) {
        if (!argv[i]) continue;
        if (std::strcmp(argv[i], "--no-file") == 0) {
            noFile = true;
        }
        else if (std::strncmp(argv[i], "--seed=", 7) == 0) {
            const char* s = argv[i] + 7;
            char* endp = nullptr;
            unsigned long v = std::strtoul(s, &endp, 10);
            if (!endp || *endp != '\0') return false;
            hasSeed = true;
            seed = static_cast<uint32_t>(v);
        }
    }

    // 3) Ask interactively if NN not provided
    if (NN == 0) {
        std::cout << "Enter NN (number of iterations, NN >= 2): ";
        if (!(std::cin >> NN)) return false;
    }

    return true;
}

// ------------------------- RNG helpers -------------------------
// rand01(): float in [0..1]
static uint32_t g_seed = 1u;

static void rng_seed(uint32_t s)
{
    g_seed = (s == 0u ? 1u : s);
}

static float rand01()
{
    // LCG (детермінований, кросплатформений)
    g_seed = 1664525u * g_seed + 1013904223u;
    return (g_seed & 0xFFFFFFu) / static_cast<float>(0xFFFFFFu);
}

// ------------------------- Small helpers -------------------------
static void push_shift(float arr[], int n, float value)
{
    // shift left, push value to the end
    for (int i = 0; i < n - 1; ++i) arr[i] = arr[i + 1];
    arr[n - 1] = value;
}

// Online mean/variance (Welford) to avoid storing NN values.
struct OnlineStats
{
    long long n = 0;
    double mean = 0.0;
    double m2 = 0.0; // sum of squares of differences from the current mean

    void add(double x)
    {
        n++;
        double delta = x - mean;
        mean += delta / n;
        double delta2 = x - mean;
        m2 += delta * delta2;
    }

    double variance_sample() const
    {
        if (n <= 1) return 0.0;
        double var = m2 / (n - 1);
        return (var < 0.0 ? 0.0 : var);
    }

    double sigma_sample() const
    {
        return std::sqrt(variance_sample());
    }
};

// ------------------------- Main -------------------------
int main(int argc, char* argv[])
{
    long NN = 0;
    bool noFile = false;
    bool hasSeed = false;
    uint32_t seed = 0;

    if (!parse_inputs(argc, argv, NN, noFile, hasSeed, seed)) {
        std::cerr << "Invalid arguments.\n"
            << "Usage: ConsoleApplication1.exe <NN> [--seed=<uint>] [--no-file]\n";
        return 1;
    }

    if (NN < 2) {
        std::cerr << "Error: NN must be >= 2.\n";
        return 1;
    }

    // Seed RNG (for reproducible tests/benchmarks)
    if (hasSeed) rng_seed(seed);
    else rng_seed(static_cast<uint32_t>(std::time(nullptr)));

    // ------------------------- Model parameters -------------------------
    // NOTE: Параметри можна змінювати у цьому блоці (Maintainability вимога).
    const float dt = 1.0f;

    const float TW = 1.0f;  // time constant for water/distillate effect
    const float TNAGR = 2.0f;  // time constant for load effect
    const float TGAZ = 1.0f;  // time constant for gas temperature effect

    const float kW = 10.0f;
    const float kNAGR = 20.0f;
    const float kGAZ = 5.0f;

    // Std-dev (sigma) of input disturbances (fluctuations)
    const float sigmaW = 1.0f;
    const float sigmaNAGR = 1.0f;
    const float sigmaGAZ = 1.0f;

    // ------------------------- Initialize buffers -------------------------
    float sdW[N]{};     // fluctuations of distillate amount (input)
    float sdNAGR[N]{};  // fluctuations of load (input)
    float sdGAZ[N]{};   // fluctuations of gas temperature (input)

    float wW[N]{}, wNAGR[N]{}, wGAZ[N]{};
    buildImpulse(TW, kW, dt, wW, N);
    buildImpulse(TNAGR, kNAGR, dt, wNAGR, N);
    buildImpulse(TGAZ, kGAZ, dt, wGAZ, N);

    // ------------------------- Simulation -------------------------
    OnlineStats stats;

    auto t0 = std::chrono::high_resolution_clock::now();

    for (long i = 0; i < NN; ++i) {
        // Generate new fluctuation samples (mean ~ 0)
        const float newW = normal12(sigmaW, rand01);
        const float newNAGR = normal12(sigmaNAGR, rand01);
        const float newGAZ = normal12(sigmaGAZ, rand01);

        push_shift(sdW, N, newW);
        push_shift(sdNAGR, N, newNAGR);
        push_shift(sdGAZ, N, newGAZ);

        // Each influence is convolution of last N samples with impulse response
        const float yW = svertka(sdW, wW, N);
        const float yNAGR = svertka(sdNAGR, wNAGR, N);
        const float yGAZ = svertka(sdGAZ, wGAZ, N);

        const float dT = yW + yNAGR + yGAZ;  // результуюче відхилення температури
        stats.add(dT);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t1 - t0).count();

    const double mean = stats.mean;
    const double sigma = stats.sigma_sample();

    // ------------------------- Output -------------------------
    std::cout << "NN = " << NN << "\n";
    std::cout << "mpar = " << mean << "\n";
    std::cout << "sig_tpar = " << sigma << "\n";
    std::cout << "time_ms = " << elapsed_ms << "\n";
    if (hasSeed) std::cout << "seed = " << seed << "\n";

    if (!noFile) {
        // Security: write only to local file name (no user-provided path)
        const char* outName = "tpar.rez";
        FILE* p1 = std::fopen(outName, "a");
        if (!p1) {
            std::cerr << "Warning: cannot open output file '" << outName << "' for writing.\n";
            return 0; // не аварійно завершуємось
        }

        std::fprintf(p1, "NN=%ld mpar=%.10f sig_tpar=%.10f time_ms=%.3f", NN, mean, sigma, elapsed_ms);
        if (hasSeed) std::fprintf(p1, " seed=%u", static_cast<unsigned>(seed));
        std::fprintf(p1, "\n");
        std::fclose(p1);
    }

    return 0;
}
