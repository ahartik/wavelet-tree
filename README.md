This is my collection of wavelet tree implementations. It is work (hopefully) towards my master's thesis.

FastBitVector
====================
Uncompressed bitvector with rank and select.
- Rank is a modification Sebastiano Vigna's rank9 from 'Broadword Implementation of Rank/Select Queries.' with less space overhead.
- Select uses its own sampling + binary search on rank-superblocks + linear search on rank-blocks.
FastBitVector Operations
-----------------------------
- FastBitVector::rank(size_t pos, bool bit)
  * Returns number of positions i < pos with vec[i] == bit
- FastBitVector::select(size_t rank, bool bit)
  * Returns smallest i with rank(i, bit) == rank

SparseBitVector
===================
Sparse bitvector 
- At the moment uses 5m + n/16 bits, in the future perhaps ~ 1.92m + m\*log2(n/m) bits like Sadakane.
- D. Okanohara, K. Sadakane: 'Practical Entropy-Compressed Rank/Select Dictionary', Proceedings of ALENEX 2007.

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
- Array of different sized balanced wavelet trees.
- Operations in O(log value)

Benchmarks
=================
On my computer (i7 2600k 4.5ghz) with popcount instruction
- FastBitVector::rank < 30ns
- FastBitVector::select < 40ns
- BalancedWaveletTree::rankLE < 1000ns
- SkewedWaveletTree::rankLE < 170ns with skewed input.

See fast-bit-vector_benchmark.cpp and wavelet_benchmark.cpp.
