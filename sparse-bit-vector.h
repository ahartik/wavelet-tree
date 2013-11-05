#pragma once

#include "fast-bit-vector.h"
#include <iterator>
#include <stdint.h>
#include <cmath>
#include <cstring>

class SparseBitVector {
 public:
  // Empty constructor
  SparseBitVector()
    : w_(4),
      pop_(0),
      size_(0),
      low_bits_(nullptr) { }

  template<typename It>
  SparseBitVector(It begin, It end) {
    init(begin, end);
  }
  SparseBitVector(const std::vector<bool>& vec) {
    std::vector<size_t> pos;
    for (size_t i = 0; i < vec.size(); ++i) {
      if (vec[i]) pos.push_back(i);
    }
    init(pos.begin(), pos.end());
  }

  SparseBitVector(SparseBitVector&& o) : SparseBitVector() {
    swap(*this, o);
  }

  const SparseBitVector& operator=(SparseBitVector&& o) {
    swap (*this, o);
    return *this;
  }

  size_t rank(size_t pos, bool bit) const {
    if (pos == 0) return 0;
    uint64_t mask = (1LL << w_) - 1;
    size_t high = pos >> w_;
    size_t low = pos & mask;
    size_t y = high_bits_.select(high, 0);
    size_t x = y - high;
    for (;high_bits_[y] == 1; x++, y++) {
      size_t l = this->low(x);
      if (l >= low) {
        // if (l == low) ++x;
        break;
      }
    }
    return bit ? x : pos - x;
  }

  // Compatibility function, for b = 0 binary search is used.
  size_t select(size_t rnk, bool b) const {
    if (rnk == 0) return 0;
    if (b == 0)
    {
      size_t left = 0;
      size_t right = size_;
      while (left != right - 1) {
        size_t c = (left + right) / 2;
        if (rank(c, 0) < rnk) {
          left = c;
        } else {
          right = c;
        }
      }
      return left + 1;
    }
    return select1(rnk);
  }
  size_t select1(size_t rank) const {
    if (rank == 0) return 0;
    return ((high_bits_.select(rank, 1) - rank) << w_) + low(rank-1) + 1;
  }

  size_t bitSize() const {
    return w_ * pop_ + high_bits_.bitSize();
  }

  size_t size() const {
    return size_;
  }

 private:
  template<typename It>
  void init(It begin, It end) {
    size_t m = std::distance(begin, end);
    pop_ = m;
    It last = end;
    size_t n = 1 + *(--last);
    size_ = n;
    // TODO calculate better w
    w_ = 4;
    size_t low_bits_size = 1 + w_ * m / 64;
    low_bits_ = new uint64_t[low_bits_size];
    memset(low_bits_, 0, 8 * low_bits_size);
    size_t j = 0;
    size_t i = 0;
    size_t offset = 0;
    uint64_t mask = (1LL << w_) - 1;
    size_t z = 1 + log2(n) - w_;
    std::vector<bool> high_bits(m + (1 << z));
    for (It it = begin; it != end; ++it) {
      uint64_t pos = *it;
      assert(pos < n);
      low_bits_[j] |= (pos & mask) << offset;
      offset += w_;
      if (offset == 64) {
        offset = 0;
        ++j;
      }
      assert(offset < 64);
      assert((pos & mask) == low(i));

      uint64_t high = pos >> w_;
      high_bits[high + i] = 1;
      ++i;
    }
    high_bits_ = FastBitVector(high_bits);
  }
  int low(size_t i) const {
    size_t p = i * w_ / 64;
    size_t o = i * w_ % 64;
    uint64_t mask = (1LL << w_) - 1;
    return (low_bits_[p] >> o) & mask;
  }
  int w_;
  size_t pop_;
  size_t size_;
  uint64_t* low_bits_;
  FastBitVector high_bits_;
 public:
  friend void swap(SparseBitVector& a, SparseBitVector& b) {
    using std::swap;
    swap(a.w_, b.w_);
    swap(a.pop_, b.pop_);
    swap(a.size_, b.size_);
    swap(a.low_bits_, b.low_bits_);
    swap(a.high_bits_, b.high_bits_);
  }
};
