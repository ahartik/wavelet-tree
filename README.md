FastBitVector
====================
This is an implementation of a simple bitvector from paper
"Fast, Small, Simple Rank/Select on Bitmaps" by Navarro and Providel.
Uses about 3% extra space for support structures.
FastBitVector Operations
-----------------------------
- FastBitVector::rank(size_t pos, bool bit)
  * Returns number of positions i < pos with vec[i] == bit
- FastBitVector::select(size_t rank, bool bit)
  * Returns smallest i with rank(i, bit) == rank


Wavelet trees
===========================
Wavelet Operations
-----------------
- \*Wavelet::rank(size_t pos, intmax_t value)
  * Returns number of positions i < pos with wt[i] == value
  * O(log n)
- \*Wavelet::rankLE(size_t pos, intmax_t value)
  * Returns number of positions i < pos with wt[i] <= value
  * O(log n)
- Some kind of select support will be added.

BalancedWavelet
-----------------
- Perfectly balanced tree for fixed number of bits per item.
- Operations in O(log n)

SkewedWavelet
-----------------
- Collection of different sized balanced wavelet trees.
- Operations in O(log value)

Benchmarks
=================
On my computer (i7 2600k 4.5ghz) with popcount instruction
- FastBitVector::rank < 30ns
- FastBitVector::select < 75ns
- BalancedWaveletTree::rankLE < 1000ns
- SkewedWaveletTree::rankLE < 170ns with skewed input.

See fast-bit-vector_benchmark.cpp and wavelet_benchmark.cpp.

TODO
=================
- Implement some kind of compressed bitvector, perhaps also from 
  "Fast, Small, Simple Rank/Select on Bitmaps"

FastBitVector
-------------------
- Get rid of std::vector<bool> in constructors - allow in-place construction:
  * FastBitVector::append OR
  * FastBitVector(PlainBitVector&& moved)
- Store small vectors in place of members (like small string optimization)

BalancedWaveletTree
-------------------
- Implement select for completeness.
- Function for listing all values <=m

SkewedWaveletTree
-------------------
- Try construction without all that extra space.
