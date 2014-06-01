#pragma once

#include <vector>

#include "fast-bit-vector.h"
#include "sdsl/rrr_vector.hpp"

class RRRBitVector {
 public:
  RRRBitVector() {
  }
  RRRBitVector(const RRRBitVector& o) = delete;

  RRRBitVector(const std::vector<bool>& vec) {
    sdsl::bit_vector bv(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
      bv[i] = vec[i];
    }
    vec_ = vector_type(bv);
    rank_.set_vector(&vec_);
    select0_.set_vector(&vec_);
    select1_.set_vector(&vec_);
  }

  RRRBitVector(RRRBitVector&& o)
      : vec_(std::move(o.vec_)),
        rank_(&vec_),
        select0_(&vec_),
        select1_(&vec_) { }
  const RRRBitVector& operator=(RRRBitVector&& o) {
    vec_ = std::move(o.vec_);
    rank_.set_vector(&vec_);
    select0_.set_vector(&vec_);
    select1_.set_vector(&vec_);
    return *this;
  }

  bool operator[](size_t i) const {
    return vec_[i];
  }
  size_t rank(size_t i, bool b) const {
    size_t r = rank_(i);
    if (!b) return i - r;
    return r;
  }
  size_t select(size_t i, bool b) {
    if (i == 0) return 0;
    if (b) return select1_(i) + 1;
    return select0_(i) + 1;
  }
  friend void swap(RRRBitVector& a, RRRBitVector& b) {
    a.vec_.swap(b.vec_);
    a.rank_.set_vector(&a.vec_);
    a.select0_.set_vector(&a.vec_);
    a.select1_.set_vector(&a.vec_);
    b.rank_.set_vector(&b.vec_);
    b.select0_.set_vector(&b.vec_);
    b.select1_.set_vector(&b.vec_);
  }
  size_t size() const {
    return vec_.size();
  }
  size_t bitSize() const {
    size_t size = 8 * (
        sdsl::size_in_bytes(vec_) + 
        sdsl::size_in_bytes(rank_) + 
        sdsl::size_in_bytes(select0_) + 
        sdsl::size_in_bytes(select1_));
    return size;
  }
 private:
  typedef sdsl::rrr_vector<63> vector_type;
  vector_type vec_;
  typename vector_type::rank_1_type rank_;
  typename vector_type::select_0_type select0_;
  typename vector_type::select_1_type select1_;
};
