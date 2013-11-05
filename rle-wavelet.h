#ifndef RLE_WAVELET_H
#define RLE_WAVELET_H

#include <algorithm>

#include "fast-bit-vector.h"
#include "sparse-bit-vector.h"
#include "skewed-wavelet.h"

class RLEWavelet {
 public:
  template<typename It>
  RLEWavelet(It begin, It end) {
    if (begin == end) return;
    typedef typename std::remove_reference<decltype(*begin)>::type ValueInt;
    size_t run = 0;
    std::vector<ValueInt> head;
    std::vector<size_t> run_end;
    head.push_back(*begin);
    size_t i = 0;
    for (It it = begin; it != end; ++it, ++i) {
      if (*it != head.back()) {
        ValueInt last = head.back();
        run_end.push_back(i);
        if (run_.size() <= last) {
          run_.resize(last + 1);
        }
        run_[last].push_back(run);
        head.push_back(*it);
        run = 0;
      }
      run++;
    }
    {
      ValueInt last = head.back();
      run_end.push_back(i);
      if (run_.size() <= last) {
        run_.resize(last + 1);
      }
      run_[last].push_back(run);
      run = 0;
    }
    assert(head.size() == run_end.size());
    run_end_ = SparseBitVector(run_end.begin(), run_end.end());

    for (size_t i = 0; i < run_.size(); ++i) {
      std::vector<size_t> nw(run_[i].size());
      size_t total = 0;
      for (size_t j = 0; j < run_[i].size(); ++j) {
        total += run_[i][j];
        nw[j] = total;
      }
      run_[i] = nw;
    }
    head_ = BalancedWavelet(head.begin(), head.end());
  }

  size_t rank(size_t pos, intmax_t value) const {
    if (pos == 0) return 0;
    size_t rpos = headPos(pos);
    BalancedWavelet::Iterator it(head_);
    bool eq = true;
    size_t hrank = rpos;
    for (;;) {
      bool bit = value >= it.splitValue();
      if (it[hrank] != bit) eq = false;
      hrank = it.rank(hrank, bit);
      if (it.isLeaf()) {
        break;
      }
      it = it.child(bit);
    }
    size_t begin = hrank == 0 ? 0 : run_[value][hrank - 1];
    size_t run_start = 0;
    if (!eq) {
      return begin;
    }
    if (rpos != 0) run_start = run_end_.select1(rpos) - 1;
    size_t end = pos - run_start;
    return begin + end;
  }

  size_t rankLE(size_t pos, intmax_t value) const {
    BalancedWavelet::Iterator it(head_);
    if (pos == 0) return 0;
    size_t rpos = headPos(pos);
    bool lt = true;
    size_t begin = rankLE(it, rpos, value, &lt);
    size_t run_start = 0;
    if (!lt) {
      return begin;
    }
    if (rpos != 0) run_start = run_end_.select1(rpos) - 1;
    size_t end = pos - run_start;
    return begin + end;
  }

  size_t size() const {
    return run_end_.size();
  }

  size_t bitSize() const {
    size_t total = head_.bitSize();
    total += run_end_.bitSize();
    std::cout << "run_end_.bitSize() = " << run_end_.bitSize() << "\n";
    for (size_t i = 0; i < run_.size(); i++) {
      total += run_[i].size() * sizeof(size_t) * 8;
    }
    return total;
  }

  int64_t operator[](size_t i) const {
    return head_[headPos(i)];
  }

 private:
  size_t headPos(size_t pos) const {
    return run_end_.rank(pos + 1, 1);
  }

  // Internal recursive function for rankLE traversal.
  // *lt is set to false if head_[pos] > value.
  size_t rankLE(BalancedWavelet::Iterator it, size_t pos, intmax_t value, bool* lt) const {
    if (it.count() == 0) return 0;
    size_t pos1 = it.rank(pos, 1);
    size_t pos0 = pos - pos1;
    intmax_t split = it.splitValue();
    bool b = value >= split;
    bool itb = it[pos];
    // All values in the future will be smaller than value
    if (itb < b) lt = nullptr;
    if (lt != nullptr && b < itb) *lt = false;
    if (it.isLeaf()) {
      size_t ret = 0;
      if (b && pos1 != 0) {
        ret += run_[split][pos1 - 1];
      }
      if (pos0 != 0) {
        ret += run_[split-1][pos0 - 1];
      }
      return ret;
    }
    if (b) {
      size_t ret = 0;
      // Only carry lt to larger child.
      ret += rankLE(it.child(0), pos0, value, nullptr);
      ret += rankLE(it.child(1), pos1, value, lt);
      return ret;
    } else {
      return rankLE(it.child(0), pos0, value, lt);
    }
  }

  SparseBitVector run_end_;
  std::vector<std::vector<size_t>> run_;
  BalancedWavelet head_;
};
#endif
