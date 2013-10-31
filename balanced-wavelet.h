#ifndef BALANCED_WAVELET_H
#define BALANCED_WAVELET_H

#include "fast-bit-vector.h"
#include <stdint.h>
#include <cassert>
#include <iostream>

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
  explicit BalancedWaveletEncoder(int b) : bits(b) { }

  BalancedWaveletEncoder(const BalancedWaveletEncoder& o) = delete;
  BalancedWaveletEncoder(BalancedWaveletEncoder&& o) {
    bits = o.bits;
    croot = o.croot;
    o.croot.vec.clear();
    o.croot.child[0] = nullptr;
    o.croot.child[1] = nullptr;
  }
  void append(intmax_t value) {
    ConstructNode* add = &croot;
    assert (value < (1LL<<bits));
    for (int i = bits-1; i > 0; --i) {
      bool b = (value >> i) & 1;
      if (add->child[b] == nullptr) {
        add->child[b] = new ConstructNode;
      }
      add->vec.push_back(b);
      add = add->child[b];
    }
    bool b = value & 1;
    add->vec.push_back(b);
  }
 private:
  ConstructNode croot;
  int bits;
  friend class BalancedWavelet;
};

class BalancedWavelet {
  struct Node {
    Node(BalancedWaveletEncoder::ConstructNode& n) : vec(n.vec) {
      for (int i = 0; i < 2; ++i) {
        if (n.child[i] != nullptr) {
          child[i] = new Node(*n.child[i]);
          n.vec.clear();
        }
      }
    }
    Node(Node&& n) : vec(std::move(n.vec)) {
      child[0] = n.child[0];
      child[1] = n.child[1];
      n.child[0] = nullptr;
      n.child[1] = nullptr;
    }
    ~Node() {
      for (int i = 0; i < 2; ++i)
        delete child[i];
    }
    FastBitVector vec;
    Node *child[2] = {nullptr, nullptr};
  };
 public:
  static const size_t npos = -1;
  BalancedWavelet(BalancedWaveletEncoder&& enc) 
      : root_(enc.croot), bits_(enc.bits)  {
  }

  template<typename It>
  BalancedWavelet(It begin, It end, int bits) 
      : BalancedWavelet(BalancedWaveletEncoder(begin, end, bits)) {
  }
  BalancedWavelet(const BalancedWavelet& o) = delete;
  BalancedWavelet(BalancedWavelet&& o) 
      : root_(std::move(o.root_)), bits_(o.bits_) { }

  size_t rank(size_t pos, int64_t value) const {
    assert(pos <= size());
    const Node* n = &root_;
    for (int i = bits_-1; i >= 0; --i) {
      if (n == nullptr) return 0;
      bool b = (value >> i) & 1;
      pos = n->vec.rank(pos, b);
      n = n->child[b];
    }
    return pos;
  }

  size_t rankLE(size_t pos, int64_t value) const {
    assert(pos <= size());
    size_t ret = 0;
    const Node* n = &root_;
    for (int i = bits_-1; i >= 0; --i) {
      if (n == nullptr) return ret;
      bool b = (value >> i) & 1;
      size_t np = n->vec.rank(pos, b);
      if (b) {
        ret += pos - np;
        assert (pos >= np);
        // std::cout << "ret = " << ret << std::endl;
      }
      pos = np;
      n = n->child[b];
    }
    return ret + pos;
  }

  size_t select(size_t rank, int64_t value) const {
    return select(&root_, rank, value, bits_-1);
  }
  
  size_t size() const {
    return root_.vec.size();
  }
  size_t bitSize() const {
    return bitSize(&root_) + 8 * sizeof(BalancedWavelet);
  }
 private:
  size_t bitSize(const Node* node) const {
    size_t ret = node->vec.extra_bits() + node->vec.size() + 8 * sizeof(Node);
    for (int i = 0; i < 2; ++i) {
      if (node->child[i] != nullptr) {
        ret += bitSize(node->child[i]);
      }
    }
    return ret;
  }

  size_t select(const Node* node, size_t rank, long value, int i) const {
    if (rank == 0) return 0;
    bool b = (value >> i) & 1;
    if (i == 0) {
      return node->vec.select(rank, b);
    }
    Node* next = node->child[b];
    if (next == nullptr) return npos;
    rank = select(next, rank, value, i - 1);
    return node->vec.select(rank, b);
  }

  Node root_;
  int bits_;
};

#endif
