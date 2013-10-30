#include "balanced-wavelet.h"
#include "skewed-wavelet.h"

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
using namespace std;

TEST(BalancedWaveletTest, Rank) {
  vector<int> v = {4,2,3,1,2,3,4,5};
  BalancedWavelet wt(v.begin(), v.end(), 3);
  EXPECT_EQ(0, wt.rank(0, 4));
  EXPECT_EQ(0, wt.rank(2, 5));

  EXPECT_EQ(1, wt.rank(1, 4));
  EXPECT_EQ(1, wt.rank(2, 2));
  EXPECT_EQ(2, wt.rank(5, 2));
}

TEST(BalancedWaveletTest, Select) {
  vector<int> v = {4,2,3,1,2,3,4,5};
  BalancedWavelet wt(v.begin(), v.end(), 3);
  EXPECT_EQ(1, wt.select(1, 4));
  EXPECT_EQ(2, wt.select(1, 2));
  EXPECT_EQ(3, wt.select(1, 3));
  EXPECT_EQ(4, wt.select(1, 1));
  EXPECT_EQ(5, wt.select(2, 2));
}
TEST(BalancedWaveletTest, RankLE) {
  vector<int> v = {4,2,3,1,2,3,4,5};
  BalancedWavelet wt(v.begin(), v.end(), 3);
  EXPECT_EQ(0, wt.rankLE(0, 4));
  EXPECT_EQ(2, wt.rankLE(2, 5));

  EXPECT_EQ(2, wt.rankLE(3, 3));
  EXPECT_EQ(3, wt.rankLE(5, 2));
  EXPECT_EQ(v.size(), wt.rankLE(v.size(), 7));
}

TEST(SkewedWaveletTest, Rank) {
  vector<int> v = {4,2,3,1,2,3,4,5};
  SkewedWavelet wt(v.begin(), v.end());
  EXPECT_EQ(0, wt.rank(0, 4));
  EXPECT_EQ(0, wt.rank(2, 5));

  EXPECT_EQ(1, wt.rank(1, 4));
  EXPECT_EQ(1, wt.rank(2, 2));
  EXPECT_EQ(2, wt.rank(5, 2));
}

TEST(SkewedWaveletTest, RankLE) {
  vector<int> v = {4,2,3,1,2,3,4,5};
  SkewedWavelet wt(v.begin(), v.end());
  EXPECT_EQ(0, wt.rankLE(0, 4));
  EXPECT_EQ(2, wt.rankLE(2, 5));

  EXPECT_EQ(2, wt.rankLE(3, 3));
  EXPECT_EQ(3, wt.rankLE(5, 2));
  EXPECT_EQ(v.size(), wt.rankLE(v.size(), 7));
}
