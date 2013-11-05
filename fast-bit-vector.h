#ifndef FASTBITVECTOR_H
#define FASTBITVECTOR_H // Small and simple bitvector.
// Designed to be fast and compact.
// Supports select and rank queries.

#include <cstddef>
#include <vector>

using std::size_t;

// Select for single word.
//  v: Input value to find position with rank r.
//  r: bit's desired rank [1-64].
// Returns: First index i so that r == rank(v,i)
int WordSelect(unsigned long v, int r);

// Plain bit-vector.
// Supports appending: use to construct FastBitVector

class PlainBitVector {
  static const int WordBits = 8 * sizeof(long);
  friend class FastBitVector;
  public:
  void push_back(bool b);
  private:
  size_t size_;
  size_t cap_;

};

class FastBitVector {
  static const unsigned RankSample = 512;
  static const unsigned SelectSample = 8192 * 4;
  static const unsigned WordBits = 8 * sizeof(long);
 public:
  // Empty constructor.
  FastBitVector();
  explicit FastBitVector(const std::vector<bool>& data);
  FastBitVector(const FastBitVector& other) = delete;
  FastBitVector(FastBitVector&& other);
  const FastBitVector& operator=(FastBitVector&& other);

  bool operator[](size_t pos) const {
    size_t i = pos / WordBits;
    int offset = pos % WordBits;
    return (bits_[i] >> offset) & 1;
  }
  // Number of positions < pos set with bit_value.
  size_t rank(size_t pos, bool bit_value) const {
    size_t block = pos / RankSample;
    size_t sum = rank_samples_[block];

    size_t word;
    for (word = block * RankSample / WordBits;
         (word + 1) * WordBits <= pos;
         ++word) {
      sum += __builtin_popcountll(bits_[word]);
    }
    size_t first_bits = pos - word * WordBits;
    // Add first bits from the last word.
    unsigned long mask = (1LL << first_bits) - 1LL;
    sum += __builtin_popcountll(bits_[word] & mask);
    if (bit_value == 0) return pos - sum;
    return sum;
  }
  // Returns smallest position pos so that rank(pos,bit) == idx
  size_t select(size_t idx, bool bit) const;

  size_t size() const {
    return size_;
  }
  size_t count(bool bit) const {
    if (bit) return popcount_;
    return size() - popcount_;
  }
  size_t extra_bits() const;
  size_t bitSize() const {
    return size() + extra_bits();
  }

  ~FastBitVector();
  friend void swap(FastBitVector& a, FastBitVector& b);
 private:

  size_t size_;
  size_t popcount_;

  unsigned long* bits_;

  size_t* rank_samples_;
  size_t* select_samples_[2];
};

#endif
