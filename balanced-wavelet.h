#ifndef BALANCED_WAVELET_H
#define BALANCED_WAVELET_H

#include "fast-bit-vector.h"
#include <stdint.h>
#include <cassert>
#include <iostream>

// TODO: use less extra space in initialization.

class BalancedWavelet {
  struct ConstructNode {
    std::vector<bool> vec;
    ConstructNode *child[2] = {nullptr, nullptr};
    ~ConstructNode() {
      delete child[0];
      delete child[1];
    }
  };
  struct Node {
    Node(ConstructNode& n) : vec(n.vec) {
      for (int i = 0; i < 2; ++i) {
        if (n.child[i] != nullptr) {
          child[i] = new Node(*n.child[i]);
        }
      }
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
  template<typename It>
  BalancedWavelet(It begin, It end, int bits) {
    std::vector<std::vector<bool>> initial;
    bits_ = bits;
    ConstructNode croot;
    for (It it = begin; it != end; ++it) {
      // consider first bits only.
      ConstructNode* add = &croot;
      int64_t value = *it;
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

    root_ = new Node(croot);
  }
  size_t rank(size_t pos, int64_t value) const {
    assert(pos <= size());
    Node* n = root_;
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
    Node* n = root_;
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
    return select(root_, rank, value, bits_-1);
  }
  
  size_t size() const {
    return root_->vec.size();
  }

  ~BalancedWavelet() {
    delete root_;
  }
 private:

  size_t select(Node* node, size_t rank, long value, int i) const {
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

  Node *root_;
  int bits_;
};

#endif
