#include "fast-bit-vector.h"

#include <cassert>
#include <cstring>
#include <cstdint>

uint8_t BytePop[256];
uint8_t ByteSelect[256][8];

int FBVInitTables() {
  for (int x = 0; x < 256; ++x) {
    BytePop[x] = __builtin_popcountll(x);
  }
  for (int x = 0; x < 256; ++x) {
    int r = 0;
    for (int i = 0; i < 8; ++i) {
      if ((x>>i)&1) {
        ByteSelect[x][r] = i+1;
        r++;
      }
    }
  }
  return 0;
}

FastBitVector::FastBitVector() {
  size_ = 0;
  bits_ = nullptr;
  rank_samples_ = nullptr;
  select_samples_[0] = nullptr;
  select_samples_[1] = nullptr;
}

FastBitVector::FastBitVector(const std::vector<bool>& data) {
  size_ = data.size();
  const size_t word_count = 1 + size_ / WordBits;
  bits_ = new unsigned long[word_count];
  memset(bits_, 0, word_count * sizeof(long));
  // Init bits.
  popcount_ = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    popcount_ += data[i];
    long p = i / WordBits;
    long o = i % WordBits;
    bits_[p] |= long(data[i]) << o;
  }

  // Init rank samples.
  rank_samples_ = new RankBlock[1 + size_ / RankSample];
  rank_samples_[0].abs = 0;
  size_t sum = 0;
  for (size_t i = 0; i <= size_/RankSample; i++) {
    uint64_t sub_block[6] = {};
    for (size_t j = 0; j < RankSample; ++j) {
      size_t idx = i * RankSample + j;
      if (idx >= data.size()) break;
      sum += data[idx];
      sub_block[j / RankSubSample] += data[idx];
    }
    for (int j = 1; j < 6; ++j) {
      sub_block[j] += sub_block[j-1];
    }
    // Put in reverse order to remove branch for sub_block = 0
    rank_samples_[i].rel = 
        (sub_block[0] << 44) + 
        (sub_block[1] << 33) + 
        (sub_block[2] << 22) + 
        (sub_block[3] << 11) + 
        (sub_block[4] << 00);
    rank_samples_[1 + i].abs = sum;
  }

  // Init select samples.
  select_samples_[1] = new uint32_t[2 + popcount_ / SelectSample];
  select_samples_[0] = new uint32_t[2 + (size_ - popcount_) / SelectSample];
  size_t sums[2] = {0, 0};
  size_t idx[2] = {1, 1};
  select_samples_[0][0] = select_samples_[1][0] = 0;
  for (size_t i = 0; i < size_; i++) {
    int bit = data[i];
    sums[bit]++;
    if (sums[bit] % SelectSample == 0) {
      size_t id = idx[bit]++;
      size_t rs = i / RankSample;
      select_samples_[bit][id] = rs;
      if (bit)
        assert(rank_samples_[rs].abs <= id * SelectSample);
#if 0
      size_t pos = i + 1;
      size_t word = pos / WordBits;
      // Set to point to the word.
      select_samples_[bit][id] = word * WordBits;
      // Correction info: how many bits are set afterwards.
      unsigned long wmask = (1L << (pos % WordBits)) - 1L;
      unsigned long counted = bit ? wmask & bits_[word] :
          (wmask & ~bits_[word]);
      select_samples_[bit][id] |= __builtin_popcountll(counted);
#endif
    }
  }
  select_samples_[0][idx[0]] = 1 + (size_ / RankSample);
  select_samples_[1][idx[1]] = 1 + (size_ / RankSample);
}

FastBitVector::FastBitVector(FastBitVector&& other) 
    : size_(0),
      popcount_(0),
      bits_(nullptr),
      rank_samples_(nullptr),
      select_samples_{nullptr,nullptr} {
  swap(*this, other);
}

const FastBitVector& FastBitVector::operator=(FastBitVector&& other) {
  swap(*this, other);
  return *this;
}

#if 0
// Number of positions < pos set with bit_value.
size_t FastBitVector::rank(size_t pos, bool bit_value) const {
  rank_count++;
  assert(pos <= size_);
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
#endif

// Returns smallest position pos so that rank(pos,bit) == idx

size_t FastBitVector::extra_bits() const {
  size_t r = 1 + popcount_ / SelectSample;
  r += 2 + (size_ - popcount_) / SelectSample;
  r += 2 + size_ / RankSample;
  return r * WordBits + sizeof(FastBitVector) * 8;
}

FastBitVector::~FastBitVector() {
  delete[] rank_samples_;
  delete[] select_samples_[0];
  delete[] select_samples_[1];
  delete[] bits_;
}
void swap(FastBitVector& a, FastBitVector& b) {
  using std::swap;
  swap(a.size_, b.size_);
  swap(a.popcount_, b.popcount_);
  swap(a.bits_, b.bits_);
  swap(a.rank_samples_, b.rank_samples_);
  swap(a.select_samples_[0], b.select_samples_[0]);
  swap(a.select_samples_[1], b.select_samples_[1]);
}
