#ifndef SKEWED_WAVELET_H
#define SKEWED_WAVELET_H

#include "fast-bit-vector.h"
#include "balanced-wavelet.h"
#include <stdint.h>
#include <vector>

class SkewedWavelet {
  const int StartSize = 2;
  const int StartBits = 1;
 public:
  template<typename It>
  SkewedWavelet(It begin, It end) {
    int64_t mx = 0;
    std::vector<std::vector<int64_t>> lv(64);
    std::vector<std::vector<bool>> pick(64);
    for (It it = begin; it != end; ++it) {
      int64_t fixed;
      int lvl = Level(*it, &fixed);
      lv[lvl].push_back(fixed);
      for (int j = 0; j < lvl; ++j)
        pick[j].push_back(0);
      pick[lvl].push_back(1);
    }
    wt_.reserve(64);
    wt_pick_.reserve(64);
    for (int lvl = 0; lvl < 64; ++lvl) {
      wt_pick_.emplace_back(pick[lvl]);
    }
    for (int lvl = 0; lvl < 64; ++lvl) {
      wt_.emplace_back(lv[lvl].begin(), lv[lvl].end(), StartBits + lvl);
    }
  }

  size_t rank(size_t pos, int64_t value) const {
    int64_t fixed;
    int lvl = Level(value, &fixed);
    for (int i = 0; i < lvl; ++i) {
      pos = wt_pick_[i].rank(pos, 0);
    }
    return wt_[lvl].rank(pos, fixed);
  }
  size_t rankLE(size_t pos, int64_t value) const {
    int64_t fixed;
    int lvl = Level(value, &fixed);
    int64_t ret = 0;
    size_t op = pos;
    for (int i = 0; i < lvl; ++i) {
      size_t np = wt_pick_[i].rank(pos, 0);
      ret += pos - np;
      pos = np;
    }
    pos = wt_pick_[lvl].rank(pos, 1);
    return ret + wt_[lvl].rankLE(pos, fixed);
  }
  size_t size() const {
    return wt_pick_[0].size();
  }
 private:
  int Level(int64_t x, int64_t* fix) const {
    int64_t size = StartSize;
    int64_t start = 0;
    int64_t end = StartSize;
    for (int i = 0; i < 64; ++i) {
      if (x < end) {
        *fix = x - start;
        return i;
      }
      size *= 2;
      start = end;
      end = start + size;
    }
    assert(false);
  }
  std::vector<FastBitVector> wt_pick_;
  std::vector<BalancedWavelet> wt_;
};

#endif
