/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include <limits.h>
#include "random_mersenne_twister.h"

using namespace feasst;

TEST(Random, MersenneTwisterran) {
  RandomMersenneTwister ran(17);
  EXPECT_EQ(5435609932269944049u, ran.int64());
  EXPECT_EQ(9787598041570588038u, ran.int64());
  EXPECT_EQ(3532934917758166329u, ran.int64());
  EXPECT_EQ(1252540491364030381u, ran.int64());
  EXPECT_EQ(14517319410219508767u, ran.int64());
  EXPECT_EQ(12107216506662538791u, ran.int64());
  EXPECT_EQ(11760184867633051564u, ran.int64());
  EXPECT_EQ(10617999215365180289u, ran.int64());
  EXPECT_EQ(720583563344608273u, ran.int64());
  EXPECT_EQ(6600496021159073464u, ran.int64());
  EXPECT_EQ(17444775666130670050u, ran.int64());
  EXPECT_EQ(1107628787702882111u, ran.int64());
  EXPECT_EQ(15938763574227812315u, ran.int64());
  ran.writeRestart("tmp/rstnr3");
  RandomMersenneTwister ran2("tmp/rstnr3");
  EXPECT_EQ(16183153811218227163u, ran.int64());
// HWH : fix restart of Mersenne Twister
//  EXPECT_EQ(16183153811218227163u, ran2.int64());
//  for (int i = 0; i < 10; ++i) {
//    EXPECT_EQ(ran.int64(), ran2.int64());
//  }
}
