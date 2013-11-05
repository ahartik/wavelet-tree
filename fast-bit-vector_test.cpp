#include<gtest/gtest.h>

#include <random>
#include "fast-bit-vector.h"

TEST(FastBitVectorTest, Index) {
  std::vector<bool> v = {false, true, true, false, true};
  FastBitVector vec(v);
  for (size_t i = 0; i < v.size(); ++i) {
    EXPECT_EQ(int(vec[i]), int(v[i]));
  }
}

TEST(FastBitVectorTest, ShortRank) {
  std::vector<bool> v = {false, true, true, false, true};
  FastBitVector vec(v);
  EXPECT_EQ(0, vec.rank(0, 1));
  EXPECT_EQ(0, vec.rank(0, 0));
  EXPECT_EQ(1, vec.rank(2, 1));
  EXPECT_EQ(1, vec.rank(2, 0));
  EXPECT_EQ(2, vec.rank(3, 1));
  EXPECT_EQ(1, vec.rank(3, 0));
}

TEST(FastBitVectorTest, LongRank) {
  int n = 1<<16;
  std::vector<bool> v;
  for (int j = 0; j < n; ++j) {
    v.push_back(bool(j%2));
  }
  FastBitVector vec(v);
  for (int j = 0; j < n; ++j) {
    ASSERT_EQ(j / 2, vec.rank(j, 1)) << j;
  }
  for (int j = 0; j < n; ++j) {
    ASSERT_EQ(1 + j / 2, vec.rank(j+1, 0)) << j;
  }
}

TEST(FastBitVectorTest, WordSelect) {
  EXPECT_EQ(1, WordSelect(0xff, 1));
  EXPECT_EQ(6, WordSelect(0xff, 6));
  EXPECT_EQ(5, WordSelect(0x5555, 3));
}

TEST(FastBitVectorTest, ShortSelect) {
  std::vector<bool> v = {false, true, true, false, true};
  FastBitVector vec(v);
  EXPECT_EQ(0, vec.select(0, 1));
  EXPECT_EQ(0, vec.select(0, 0));

  EXPECT_EQ(2, vec.select(1, 1));
  EXPECT_EQ(1, vec.select(1, 0));

  EXPECT_EQ(3, vec.select(2, 1));
  EXPECT_EQ(4, vec.select(2, 0));
}

TEST(FastBitVectorTest, LongSelect) {
  int n = 1<<16;
  int m = n / 4;
  std::vector<bool> v;
  for (int j = 0; j < n; ++j) {
    v.push_back(bool(j%2));
  }
  FastBitVector vec(v);
  for (int j = 0; j < m; ++j) {
    ASSERT_EQ(j*2, vec.select(j,1)) << j;
  }
}

TEST(FastBitVectorTest, RandomRank) {
  std::mt19937_64 mt(0);
  int n = 1<<20;
  int m = n / 4;
  std::vector<bool> v;
  for (int j = 0; j < n; ++j) {
    v.push_back(mt()%2);
  }
  FastBitVector vec(v);
  int rank = 0;
  for (int j = 0; j < m; ++j) {
    ASSERT_EQ(rank, vec.rank(j, 1));
    rank += v[j];
  }
}

TEST(FastBitVectorTest, RandomSelect) {
  std::mt19937_64 mt(0);
  int n = 1<<20;
  int m = n / 4;
  std::vector<bool> v;
  for (int j = 0; j < n; ++j) {
    v.push_back(mt()%2);
  }
  FastBitVector vec(v);
  int rank[2] = {0,0};
  for (int j = 0; j < m; ++j) {
    int b = v[j];
    rank[b]++;
    ASSERT_EQ(j+1, vec.select(rank[b],b)) << j;
  }
}

TEST(FastBitVectorTest, RandomRankSelect) {
  std::mt19937_64 mt(0);
  const int rarity = 1<<10;
  int n = 1<<20;
  int m = n / 4;
  std::vector<bool> v;
  for (int j = 0; j < n; ++j) {
    v.push_back(mt()%rarity == 0);
  }
  FastBitVector vec(v);
  for (int j = 1; j <= vec.count(1); ++j) {
    int b = mt() % 2;
    int p = vec.select(j, b);
    ASSERT_EQ(vec.rank(p, b), j);
  }
}

TEST(FastBitVectorTest, Edges) {
  int n = 1<<16;
  std::vector<bool> v(n, true);
  FastBitVector vec(v);
  // This works even without zeros.
  ASSERT_EQ(0, vec.select(0, 0)); 
  // especially 0 == vec.rank(0,1)
  // and vec.rank(n,1) = n
  for (int j = 0; j <= n; ++j) {
    ASSERT_EQ(j, vec.rank(j,1)) << j;
    ASSERT_EQ(0, vec.rank(j,0)) << j;
    ASSERT_EQ(j, vec.select(j,1)) << j;
  }
}
