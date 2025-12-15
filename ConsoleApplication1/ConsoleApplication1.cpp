// rodion1_updated.cpp
// Console application for statistical modeling of steam temperature fluctuations
// (based on the original rodion1.cpp, with improved I/O, validation, and timing)
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>
#include <cstdint>

#define N 60

static float svertka(const float x[], const float w[], int n);

// Parse CLI args: <NN> [--no-file] [--seed=<uint>]
// If NN is not provided, reads from stdin.
static bool parse_inputs(int argc, char* argv[], long& NN, bool& noFile, bool& hasSeed, uint32_t& seed)
{
    NN = -1;
    noFile = false;
    hasSeed = false;
    seed = 0;

    // Scan flags
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--no-file") {
            noFile = true;
            continue;
        }
        const std::string pref = "--seed=";
        if (a.rfind(pref, 0) == 0) {
            std::string v = a.substr(pref.size());
            try {
                unsigned long tmp = std::stoul(v);
                seed = static_cast<uint32_t>(tmp);
                hasSeed = true;
            }
            catch (...) {
                return false;
            }
            continue;
        }

        // First non-flag token treated as NN (if not yet set)
        if (NN < 0) {
            try {
                size_t pos = 0;
                long val = std::stol(a, &pos);
                if (pos != a.size()) return false;
                NN = val;
            }
            catch (...) {
                return false;
            }
        }
        else {
            // Unexpected token
            return false;
        }
    }

    if (NN >= 0) return true;

    // Fallback: stdin
    std::cout << "Input NN: ";
    if (!(std::cin >> NN)) return false;
    return true;
}

int main(int argc, char* argv[])
{
    using clock = std::chrono::high_resolution_clock;
    const auto t_program_start = clock::now();

    // ----------------- Model parameters (can be changed before build) -----------------
    const float TW = 1.0f, TNAGR = 2.0f, TGAZ = 1.0f;
    const float dt = 1.0f;
    const float kW = 10.0f, kNAGR = 20.0f, kGAZ = 5.0f;

    // Standard deviations of input fluctuations
    const float sigma_W = 0.05f, sigma_NAGR = 0.05f, sigma_GAZ = 5.0f;
    // -------------------------------------------------------------------------------

    // ----------------- Input -----------------
    long NN = 0;
    bool noFile = false;
    bool hasSeed = false;
    uint32_t seed = 0;

    if (!parse_inputs(argc, argv, NN, noFile, hasSeed, seed)) {
        std::cout << "ERROR: Invalid arguments. Usage: <program> <NN> [--no-file] [--seed=<uint>]\n";
        std::cout << "ERROR: NN must be an integer.\n";
        return 1;
    }

    if (NN < 2 || NN > 1000000) {
        std::cout << "ERROR: NN must be in range [2..1000000].\n";
        return 1;
    }

    if (hasSeed) {
        std::srand(seed);
    }
    else {
        // Non-deterministic seed for typical runs
        std::srand(static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    }
    // ----------------------------------------

    // Startup time: till after argument parsing & validation
    const auto t_startup_done = clock::now();
    const auto startup_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_startup_done - t_program_start).count();

    // ----------------- Prepare I/O -----------------
    std::FILE* p1 = nullptr;
    if (!noFile) {
        p1 = std::fopen("tpar.rez", "a");
        if (!p1) {
            std::cout << "ERROR: cannot open output file tpar.rez.\n";
            return 1;
        }
    }
    // ----------------------------------------------

    float sdNAGR[N]{};
    float sdGAZ[N]{};
    float sdW[N]{};

    float wW[N]{}, wNAGR[N]{}, wGAZ[N]{};

    // Impulse response functions for 1st order links
    for (int i = 0; i < N; ++i) wW[i] = (kW / TW) * std::exp(-i * dt / TW);
    for (int i = 0; i < N; ++i) wNAGR[i] = (kNAGR / TNAGR) * std::exp(-i * dt / TNAGR);
    for (int i = 0; i < N; ++i) wGAZ[i] = (kGAZ / TGAZ) * std::exp(-i * dt / TGAZ);

    // ----------------- Simulation -----------------
    const auto t_sim_start = clock::now();

    double stpar = 0.0;   // sum(dtpar)
    double stpar2 = 0.0;  // sum(dtpar^2)

    for (long j = 0; j < NN - 1; ++j) {
        // shift histories left (older -> index 0)
        for (int i = 0; i < N - 1; ++i) {
            sdGAZ[i] = sdGAZ[i + 1];
            sdNAGR[i] = sdNAGR[i + 1];
            sdW[i] = sdW[i + 1];
        }

        // Normal-like random by sum of 12 uniform(0..1) - 6
        float sl = 0.0f;
        for (int i = 0; i < 12; ++i) sl += static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        sdGAZ[N - 1] = sigma_GAZ * (sl - 6.0f);

        sl = 0.0f;
        for (int i = 0; i < 12; ++i) sl += static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        sdNAGR[N - 1] = sigma_NAGR * (sl - 6.0f);

        sl = 0.0f;
        for (int i = 0; i < 12; ++i) sl += static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        sdW[N - 1] = sigma_W * (sl - 6.0f);

        // Output temperature fluctuation: sum of three convolutions
        float dtpar = svertka(sdGAZ, wGAZ, N) + svertka(sdNAGR, wNAGR, N) + svertka(sdW, wW, N);

        stpar += dtpar;
        stpar2 += static_cast<double>(dtpar) * static_cast<double>(dtpar);
    }

    // Statistics
    const double mtpar = stpar / static_cast<double>(NN);
    const double Dtpar = (stpar2 - stpar * stpar / static_cast<double>(NN)) / static_cast<double>(NN - 1);
    double sig_tpar = 0.0;
    bool ok = true;
    if (Dtpar >= 0.0) {
        sig_tpar = std::sqrt(Dtpar);
    }
    else {
        ok = false;
    }

    const auto t_sim_end = clock::now();
    const auto sim_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_sim_end - t_sim_start).count();
    const double us_per_iter = (static_cast<double>(sim_ms) * 1000.0) / static_cast<double>(NN);
    // ---------------------------------------------

    // ----------------- Output -----------------
    if (!ok) {
        std::cout << "ERROR: Variance became negative (Dtpar < 0). Check NN and numeric stability.\n";
        if (p1) std::fclose(p1);
        return 1;
    }

    std::cout << "Startup time (ms): " << startup_ms << "\n";
    std::cout << "Total simulation time (ms): " << sim_ms << "\n";
    std::cout << "Avg latency per iteration (us): " << us_per_iter << "\n";
    std::cout << "mpar=" << mtpar << " sig_tpar=" << sig_tpar << "\n";
    std::cout << "Finish!\n";

    if (p1) {
        std::fprintf(p1,
            "NN=%ld; mpar=%f; sig_tpar=%f; startup_ms=%lld; sim_ms=%lld; us_per_iter=%.3f; seed=%u%s\n",
            NN,
            static_cast<float>(mtpar),
            static_cast<float>(sig_tpar),
            static_cast<long long>(startup_ms),
            static_cast<long long>(sim_ms),
            us_per_iter,
            hasSeed ? seed : 0u,
            hasSeed ? "" : " (auto)"
        );
        std::fclose(p1);
    }

    return 0;
}

static float svertka(const float x[], const float w[], int n)
{
    float s = 0.0f;
    for (int i = 0; i < n; ++i) {
        s += w[i] * x[n - 1 - i];
    }
    return s;
}
