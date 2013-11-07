#include "skewed-wavelet.h"
#include "balanced-wavelet.h"
#include "rle-wavelet.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace std;

template<typename Wt>
void RankLE(int iters, int m, const char* name) {
  const size_t size = 1<<20;
  const size_t max = 1<<20;
  std::cout << name << "::rankLE(" << m << "):\n";
  using namespace std::chrono;
  std::mt19937_64 mt(0);
  std::vector<uint64_t> v;
  for (size_t j = 0; j < size; ++j) {
    v.push_back(mt()%max);
  }
  Wt wt(v.begin(), v.end());
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
  RankLE<BalancedWavelet<>>(100000, 32, "BalancedWavelet");
  RankLE<BalancedWavelet<>>(100000, 1<<10, "BalancedWavelet");
  cout << endl;
  RankLE<SkewedWavelet<>>(100000, 32, "SkewedWavelet");
  RankLE<SkewedWavelet<>>(100000, 1<<10, "SkewedWavelet");
  cout << endl;
  RankLE<RLEWavelet<BalancedWavelet<>>>(100000, 32, "RLEWavelet<BalancedWavelet>");
  RankLE<RLEWavelet<BalancedWavelet<>>>(100000, 1<<10, "RLEWavelet<BalancedWavelet>");
  cout << endl;
  RankLE<RLEWavelet<SkewedWavelet<>>>(100000, 32, "RLEWavelet<SkewedWavelet>");
  RankLE<RLEWavelet<SkewedWavelet<>>>(100000, 1<<10, "RLEWavelet<SkewedWavelet>");
}
