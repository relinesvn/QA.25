#pragma once
#include <cstdint>

constexpr int N = 60;

// Convolution (your svertka)
float svertka(const float x[], const float w[], int n);

// Generate normal-like random value using "sum of 12 uniforms - 6", scaled by sigma.
// rand01 must return float in [0..1].
float normal12(float sigma, float (*rand01)());

// Build impulse response for 1st order link: w[i] = (k/T)*exp(-i*dt/T), i=0..n-1
void buildImpulse(float T, float k, float dt, float outW[], int n);

// Compute mean and sample std-dev (sigma) for array values[0..count-1]
void computeMeanSigma(const float* values, int count, double& mean, double& sigma);
