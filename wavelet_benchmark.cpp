#include "skewed-wavelet.h"
#include "balanced-wavelet.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace std;
void SkewedRankLE(int iters, int m) {
  const size_t size = 1<<20;
  const size_t max = 1<<30;
  std::cout << "SkewedWavelet::RankLE(" << m << "):\n";
  using namespace std::chrono;
  std::mt19937_64 mt(0);
  std::binomial_distribution<int> gen(max, 4.0/max);
  std::vector<int> v;
  for (size_t j = 0; j < size; ++j) {
    v.push_back(gen(mt));
  }
  SkewedWavelet wt(v.begin(), v.end());
  std::chrono::high_resolution_clock clock;
  auto start = clock.now();
  unsigned long long total = 0;
  for (int j = 0; j < iters; ++j) {
    total = total * 178923 + 987341;
    total += wt.rankLE(total % size, m);
  }
  // Make sure compiler is not too smart.
  std::cout << "(" << total << ")\n";
  
  auto end = clock.now();
  std::cout << duration_cast<nanoseconds>(end-start).count()/iters << "ns/rank\n";
}

void BalancedRankLE(int iters, int m) {
  const size_t size = 1<<20;
  const size_t max = 1<<30;
  std::cout << "BalancedWavelet::RankLE(" << m << "):\n";
  using namespace std::chrono;
  std::mt19937_64 mt(0);
  std::binomial_distribution<int> gen(max, 4.0/max);
  BalancedWaveletEncoder enc(31);
  for (size_t j = 0; j < size; ++j) {
    enc.append(gen(mt));
  }
  BalancedWavelet wt(std::move(enc));
  std::chrono::high_resolution_clock clock;
  auto start = clock.now();
  unsigned long long total = 0;
  for (int j = 0; j < iters; ++j) {
    total = total * 178923 + 987341;
    total += wt.rankLE(total % size, m);
  }
  // Make sure compiler is not too smart.
  std::cout << "(" << total << ")\n";
  auto end = clock.now();
  std::cout << duration_cast<nanoseconds>(end-start).count()/iters << "ns/rank\n";
}

int main() {
  for (int i = 1; i < 16; ++i) {
    SkewedRankLE(100000, i);
  }
  SkewedRankLE(100000, 1<<20);
  for (int i = 1; i < 16; ++i) {
    BalancedRankLE(100000, i);
  }
}
