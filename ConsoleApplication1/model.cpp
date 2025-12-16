#include "model.h"
#include <cmath>


float svertka(const float x[], const float w[], int n)
{
    float s = 0.0f;
    for (int i = 0; i < n; ++i) {
        s += w[i] * x[n - 1 - i];
    }
    return s;
}

float normal12(float sigma, float (*rand01)())
{
    if (sigma == 0.0f) return 0.0f;

    float sl = 0.0f;
    for (int i = 0; i < 12; ++i) sl += rand01();     // sum U[0..1]
    return sigma * (sl - 6.0f);                      // ~N(0, sigma^2) approx
}

void buildImpulse(float T, float k, float dt, float outW[], int n)
{
    // minimal validation (so tests can expect deterministic behavior)
    if (n <= 0) return;

    for (int i = 0; i < n; ++i) {
        outW[i] = (k / T) * std::exp(-i * dt / T);
    }
}

void computeMeanSigma(const float* values, int count, double& mean, double& sigma)
{
    mean = 0.0;
    sigma = 0.0;
    if (!values || count <= 0) return;
    if (count == 1) { mean = values[0]; sigma = 0.0; return; }

    double sum = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < count; ++i) {
        sum += values[i];
        sum2 += static_cast<double>(values[i]) * static_cast<double>(values[i]);
    }

    mean = sum / count;
    double var = (sum2 - (sum * sum) / count) / (count - 1);  // sample variance
    if (var < 0.0) var = 0.0;                                 // clamp numeric noise
    sigma = std::sqrt(var);
}
