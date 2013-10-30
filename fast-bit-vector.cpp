#include "fast-bit-vector.h"

#include <cassert>
#include <cstring>
#include <cstdint>

static uint8_t BytePop[256];
static uint8_t ByteSelect[256][8];

static int InitTables() {
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

// Select for single word.
//  v: Input value to find position with rank r.
//  r: bit's desired rank [1-64].
// Returns: First index i so that r == rank(v,i)
int WordSelect(unsigned long v, int r) {
  static int init = InitTables();
  (void) init;
  for (size_t b = 0; b < sizeof(long); ++b) {
    int c = BytePop[v&255];
    if (c >= r) {
      return b * 8 + ByteSelect[v&255][r-1];
    }
    r -= c;
    v >>= 8;
  }
  assert(false);
  return -1;
}

static const int WordBits = 8 * sizeof(long);
static const int RankSample = 1024 * 2;
static const int SelectSample = 8192 * 4;

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
  rank_samples_ = new size_t[1 + size_ / RankSample];
  rank_samples_[0] = 0;
  size_t sum = 0;
  for (size_t i = 0; i < size_/RankSample; i++) {
    for (int j = 0; j < RankSample; ++j) {
      size_t idx = i * RankSample + j;
      sum += data[idx];
    }
    rank_samples_[1 + i] = sum;
  }

  // Init select samples.
  select_samples_[1] = new size_t[2 + popcount_ / SelectSample];
  select_samples_[0] = new size_t[2 + (size_ - popcount_) / SelectSample];
  size_t sums[2] = {0, 0};
  size_t idx[2] = {1, 1};
  select_samples_[0][0] = select_samples_[1][0] = 0;
  for (size_t i = 0; i < size_; i++) {
    int bit = data[i];
    sums[bit]++;
    if (sums[bit] % SelectSample == 0) {
      size_t id = idx[bit]++;
      size_t pos = i + 1;
      size_t word = pos / WordBits;
      // Set to point to the word.
      select_samples_[bit][id] = word * WordBits;
      // Correction info: how many bits are set afterwards.
      unsigned long wmask = (1L << (pos % WordBits)) - 1L;
      unsigned long counted = bit ? wmask & bits_[word] :
          (wmask & ~bits_[word]);
      select_samples_[bit][id] |= __builtin_popcountll(counted);
    }
  }
  select_samples_[0][idx[0]] = (size_ / WordBits) * WordBits;
  select_samples_[1][idx[1]] = (size_ / WordBits) * WordBits;
}

FastBitVector::FastBitVector(FastBitVector&& other) 
    : size_(0),
      popcount_(0),
      bits_(nullptr),
      rank_samples_(nullptr),
      select_samples_{nullptr,nullptr} {
  swap(*this, other);
}

// Number of positions < pos set with bit_value.
size_t FastBitVector::rank(size_t pos, bool bit_value) const {
  assert(pos <= size_);
  size_t block = pos / RankSample;
  size_t sum = rank_samples_[block];

  size_t word;
  for (word = block * RankSample / WordBits;
       (word + 1) * WordBits <= pos;
       ++word) {
    sum += __builtin_popcountll(bits_[word]);
  }
  long first_bits = pos - word * WordBits;
  // Add first bits from the last word.
  unsigned long mask = (1LL<<first_bits) - 1LL;
  sum += __builtin_popcountll(bits_[word] & mask);
  if (bit_value == 0) return pos - sum;
  return sum;
}

// Returns smallest position pos so that rank(pos,bit) == idx
size_t FastBitVector::select(size_t idx, bool bit) const {
  assert(idx <= count(bit));
  // flip this to operate without select samples.
#if 1
  // Start from sampling.
  size_t block = idx / SelectSample;
  size_t sample_value = select_samples_[bit][block];
  // Scan rank blocks with binary search.
  size_t word_start = sample_value / WordBits;
  size_t rank_correct = sample_value % WordBits;
  size_t word_rank = block * SelectSample - rank_correct;
  
  size_t left = (word_start * WordBits + RankSample - 1) / RankSample;

  // Scan rank blocks with binary search.
  // Look for right end using select_samples_.
  size_t sample_next = select_samples_[bit][block+1];
  size_t word_end = sample_next / WordBits;
  size_t right = (word_end * WordBits + RankSample - 1) / RankSample;
#else 
  size_t left = 0;
  size_t right = 1 + size_ / RankSample;
  size_t word_start = 0;
  size_t word_rank = 0;
#endif
  while (left + 1 < right) {
    size_t c = (left + right) / 2;
    size_t r = bit ? rank_samples_[c] : c * RankSample - rank_samples_[c];
    if (r >= idx) {
      right = c;
    } else {
      left = c;
      word_start = c * RankSample / WordBits;
      word_rank = r;
    }
  }

  // Scan words
  const size_t total_words = (size_ + WordBits - 1) / WordBits;
  for (size_t w = word_start; w < total_words; ++w) {
    size_t pop = __builtin_popcountll(bits_[w]);
    size_t r = word_rank + (bit ? pop : WordBits - pop);
    if (r >= idx) break;
    word_rank = r;
    word_start = w + 1;
  }

  // Scan bits
  size_t id = idx - word_rank;
  assert(id >= 0 && id <= WordBits);
  if (id == 0) return word_start * WordBits;
  size_t w = bits_[word_start];
  if (!bit) w = ~w;
  return word_start * WordBits + WordSelect(w, id);
}

size_t FastBitVector::extra_bits() const {
  size_t r = 1 + popcount_ / SelectSample;
  r += 2 + (size_ - popcount_) / SelectSample;
  r += 2 + size_ / RankSample;
  return r * WordBits;
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
