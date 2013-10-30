#include "skewed-wavelet.h"
#include "balanced-wavelet.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace std;
void RankLE(int iters, int m) {
  const size_t size = 1<<20;
  std::cout << "RankLE(" << m << "):\n";
  using namespace std::chrono;
  std::mt19937_64 mt(time(0));
  std::binomial_distribution<int> gen(128, 0.7);
  std::vector<int> v;
  for (size_t j = 0; j < size; ++j) {
    v.push_back(gen(mt));
  }
  SkewedWavelet wt(v.begin(), v.end());
  std::cout << "constructed\n";
  std::chrono::high_resolution_clock clock;
  auto start = clock.now();
  unsigned long long total = 0;
  for (int j = 0; j < iters; ++j) {
    total = total * 178923 + 987341;
    total += wt.rankLE(total % size, m);
  }
  std::cout << "total = " << total << endl;
  auto end = clock.now();
  long ms = duration_cast<std::chrono::milliseconds>(end-start).count();
  std::cout << ms << "ms\n";
  std::cout << duration_cast<nanoseconds>(end-start).count()/iters << "ns/rank\n";
}

int main() {
  for (int i = 1; i < 16; ++i) {
    RankLE(100000, i);
  }
}
