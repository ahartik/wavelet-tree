#ifndef BALANCED_WAVELET_H
#define BALANCED_WAVELET_H

#include "fast-bit-vector.h"
#include <stdint.h>
#include <cassert>
#include <iostream>
#include <cmath>

class BalancedWaveletEncoder {
  struct ConstructNode {
    std::vector<bool> vec;
    ConstructNode *child[2] = {nullptr, nullptr};
    ~ConstructNode() {
      delete child[0];
      delete child[1];
    }
  };
 public:
  template <typename It>
  BalancedWaveletEncoder(It begin, It end, int b) : bits(b) {
    for (It it = begin; it != end; ++it) {
      append(*it);
    }
  }

  template <typename It>
  BalancedWaveletEncoder(It begin, It end) {
    intmax_t max = 0;
    for (It it = begin; it != end; ++it) {
      if (*it > max) {
        max = *it;
      }
    }
    bits = 1 + log2(max);
    for (It it = begin; it != end; ++it) {
      append(*it);
    }
  }

  explicit BalancedWaveletEncoder(int b) : bits(b) { }

  void append(intmax_t value) {
    assert (value < (1LL<<bits));
    ConstructNode* add = &croot;
    int j = 0;
    for (int i = bits - 1; i >= 0; --i) {
      bool b = (value >> i) & 1;
      if (add->child[b] == nullptr) {
        add->child[b] = new ConstructNode;
      }
      add->vec.push_back(b);
      add = add->child[b];
      j++;
    }
  }

 private:
  friend class BalancedWavelet;
  int bits;
  ConstructNode croot;
};

class BalancedWavelet {
 public:
  static const size_t npos = -1;
  BalancedWavelet(BalancedWaveletEncoder&& enc) 
      : size_(enc.croot.vec.size()),
        bits_(enc.bits) {
    std::vector<bool> init;
    init.reserve(bits_ * enc.croot.vec.size());

    std::vector<BalancedWaveletEncoder::ConstructNode*> q;
    q.push_back(&enc.croot);
    for (size_t i = 0; i < q.size(); ++i) {
      auto* n = q[i];
      init.insert(init.end(), n->vec.begin(), n->vec.end());
      // Clear the vector to free memory.
      n->vec = std::vector<bool>();
      for (int j = 0; j < 2; ++j)
        if (n->child[j] != nullptr)
          q.push_back(n->child[j]);
    }
    delete enc.croot.child[0];
    delete enc.croot.child[1];
    enc.croot.child[0] = nullptr;
    enc.croot.child[1] = nullptr;
    FastBitVector tree(init);
    swap(tree_, tree);
  }

  template<typename It>
  BalancedWavelet(It begin, It end, int bits) 
      : BalancedWavelet(BalancedWaveletEncoder(begin, end, bits)) {
  }

  template<typename It>
  BalancedWavelet(It begin, It end) 
      : BalancedWavelet(BalancedWaveletEncoder(begin, end)) {
  }

  BalancedWavelet(const BalancedWavelet& o) = delete;

  BalancedWavelet() {
    size_ = 0;
    bits_ = 0;
  }
  BalancedWavelet(BalancedWavelet&& o) 
      : tree_(std::move(o.tree_)), size_(o.size_), bits_(o.bits_) { }

  const BalancedWavelet& operator=(BalancedWavelet&& o) {
    tree_ = std::move(o.tree_);
    size_ = o.size_;
    bits_ = o.bits_;
    o.size_ = 0;
    o.bits_ = 0;
    return *this;
  }

  class Iterator {
   public:
    Iterator(const BalancedWavelet& wt)
        : high_bits(0),
          len(wt.size_),
          offset(0),
          bit(wt.bits_ - 1),
          begin_rank(0),
          end_rank(wt.tree_.rank(wt.size_, 1)),
          level_skip(wt.size_),
          vec(&wt.tree_)
    {}

    intmax_t splitValue() const {
      return high_bits + (1LL << bit);
    }

    bool isLeaf() const {
      return bit == 0;
    }

    Iterator child(bool right) {
      return Iterator(*this, right);
    }

    bool operator[](size_t i) const {
      return (*vec)[offset + i];
    }

    size_t rank(size_t pos, bool bit) const {
      assert(pos <= len);
      size_t orank = bit ? begin_rank : offset - begin_rank;
      return vec->rank(offset + pos, bit) - orank;
    }
    size_t select(size_t idx, bool bit) const {
      size_t orank = bit ? begin_rank : offset - begin_rank;
      return vec->select(idx + orank, bit) - offset;
    }
    size_t count() const {
      return len;
    }
   private:
    Iterator(const Iterator& parent, bool right) : vec(parent.vec) {
      bit = parent.bit - 1;
      level_skip = parent.level_skip;
      if (right) {
        offset = parent.offset + level_skip + (parent.len - parent.end_rank);
        len = parent.end_rank;
        high_bits = parent.high_bits + (1LL << parent.bit);
      } else {
        offset = parent.offset + level_skip;
        high_bits = parent.high_bits;
        len = parent.len - parent.end_rank;
      }
      begin_rank = vec->rank(offset, 1);
      end_rank = vec->rank(offset + len, 1) - begin_rank;
    }
    intmax_t high_bits;
    size_t len;
    size_t offset;
    size_t bit;
    size_t begin_rank;
    size_t end_rank;
    size_t level_skip;
    const FastBitVector* vec;
    friend class BalancedWavelet;
  };

  size_t rank(size_t pos, int64_t value) const {
    assert(pos <= size());
    Iterator it(*this);
    for (;;) {
      bool bit = value >= it.splitValue();
      pos = it.rank(pos, bit);
      if (it.isLeaf()) break;
      it = it.child(bit);
    }
    return pos;
  }

  size_t rankLE(size_t pos, int64_t value) const {
    assert(pos <= size());
    size_t ret = 0;
    Iterator it(*this);
    for (;;) {
      bool bit = value >= it.splitValue();
      size_t np = it.rank(pos, bit);
      if (bit) {
        ret += pos - np;
      }
      pos = np;
      if (it.isLeaf()) break;
      else it = it.child(bit);
    }
    return pos + ret;
  }

  intmax_t operator[](size_t i) const {
    Iterator it(*this);
    while (!it.isLeaf()) {
      bool b = it[i];
      Iterator nit = it.child(b);
      i = it.rank(i, b);
      it = nit;
    }
    return it.high_bits + it[i];
  }

  size_t select(size_t rank, int64_t value) const {
    return select(Iterator(*this), rank, value);
  }

  size_t select(Iterator it, size_t rank, long value) const {
    if (it.count() == 0) return 0;
    if (rank == 0) return 0;
    bool b = value >= it.splitValue();
    if (it.isLeaf()) {
      return it.select(rank, b);
    }
    rank = select(it.child(b), rank, value);
    return it.select(rank, b);
  }

  size_t size() const {
    return size_;
  }
  size_t bitSize() const {
    return tree_.size() + tree_.extra_bits() + sizeof(*this) * 8;
  }
 private:

  FastBitVector tree_;
  size_t size_;
  int bits_;
};

#endif
