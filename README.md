This is my garden of wavelet trees. It is (hopefully) work towards my master's thesis.

FastBitVector
====================
Uncompressed bitvector with rank and select.
- Rank is a modification Sebastiano Vigna's rank9 from 'Broadword Implementation of Rank/Select Queries.' with less space overhead.
- Select uses its own sampling + binary search on rank-superblocks + linear search on rank-blocks.

SparseBitVector
===================
Sparse bitvector 
- At the moment uses 5m + n/16 bits, in the future perhaps ~ 1.92m + m\*log2(n/m) bits like Sadakane.
- D. Okanohara, K. Sadakane: 'Practical Entropy-Compressed Rank/Select Dictionary', Proceedings of ALENEX 2007.

Wavelet trees
===========================

BalancedWavelet
-----------------
- Perfectly balanced tree for fixed number of bits per item.
- Operations in O(log n)

SkewedWavelet
-----------------
- Array of different sized balanced wavelet trees.
- Operations in O(log value)

RLEWavelet<Wavelet>
-----------------
- Run length compressed Wavelet tree, uses SparseBitVector for run lengths.
- rank = O(Wavelet::rank)
- rankLE = O(Wavelet::rank + value);

Benchmarks
=================
On my computer (i7 2600k 4.5ghz) with popcnt instruction
- FastBitVector::rank < 30ns
- FastBitVector::select < 40ns
- See fast-bit-vector_benchmark.cpp and wavelet_benchmark.cpp.
