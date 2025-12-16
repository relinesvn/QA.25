#include "pch.h"
#include "CppUnitTest.h"
#include "../ConsoleApplication1/model.h"  
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace ConsoleApplication1Tests
{
    // deterministic "rand01" for tests
    static float fakeRand01()
    {
        // simple LCG-based deterministic sequence in [0..1]
        static uint32_t s = 1u;
        s = 1664525u * s + 1013904223u;
        return (s & 0xFFFFFF) / static_cast<float>(0xFFFFFF);
    }

    TEST_CLASS(UnitCoreTests)
    {
    public:

        TEST_METHOD(Svertka_KnownVector_AllOnes)
        {
            // REQ-SV-1
            float x[3] = { 1,2,3 };
            float w[3] = { 1,1,1 };
            float y = svertka(x, w, 3);
            Assert::AreEqual(6.0f, y, 1e-6f);
        }

        TEST_METHOD(Svertka_ZeroInput_ReturnsZero)
        {
            // REQ-SV-2
            float x[4] = { 0,0,0,0 };
            float w[4] = { 5,4,3,2 };
            float y = svertka(x, w, 4);
            Assert::AreEqual(0.0f, y, 1e-6f);
        }

        TEST_METHOD(Normal12_SigmaZero_ReturnsZeroAlways)
        {
            // REQ-RND-1
            for (int i = 0; i < 100; ++i) {
                float v = normal12(0.0f, fakeRand01);
                Assert::AreEqual(0.0f, v, 0.0f);
            }
        }

        TEST_METHOD(BuildImpulse_FirstElement_Equals_k_div_T)
        {
            // REQ-IMP-1
            float w[N]{};
            const float T = 2.0f, k = 10.0f, dt = 1.0f;
            buildImpulse(T, k, dt, w, N);
            Assert::AreEqual(k / T, w[0], 1e-6f);
        }

        TEST_METHOD(ComputeStats_ConstantSeries_SigmaZero)
        {
            // REQ-STAT-1
            float v[10];
            for (int i = 0; i < 10; ++i) v[i] = 7.0f;

            double mean = 0.0, sigma = 0.0;
            computeMeanSigma(v, 10, mean, sigma);

            Assert::AreEqual(7.0, mean, 1e-12);
            Assert::AreEqual(0.0, sigma, 1e-12);
        }
    };
}
