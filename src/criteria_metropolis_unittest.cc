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
#include "criteria.h"
#include "criteria_metropolis.h"

TEST(Criteria, Metropolis) {
  const double beta = 1e-15, activ = 1;
  feasst::CriteriaMetropolis c(beta, activ);
  EXPECT_EQ(1, c.accept(0., 0., "move", 0));
}

TEST(Criteria, args) {
  const double beta = 0.7;
  feasst::CriteriaMetropolis c(beta, {{"activ", "0.2345"}});
  EXPECT_NEAR(c.activ(0), 0.2345, feasst::DTOL);

  try {
    feasst::CriteriaMetropolis c(beta, {{"/not/an/arg", "error"}});
    CATCH_PHRASE("is not recognized");
  }

  auto criteria = feasst::makeCriteriaMetropolis(
    {{"beta", "0.69687"},
     {"activ", "456.3486"}});
  EXPECT_NEAR(criteria->beta(), 0.69687, feasst::DTOL);
  EXPECT_NEAR(criteria->activ(), 456.3486, feasst::DTOL);

  try {
    auto criteria = feasst::makeCriteriaMetropolis({{"not/an/arg/", "error"}});
    CATCH_PHRASE("key(beta) is required for args");
  }

  try {
    auto criteria = feasst::makeCriteriaMetropolis(
      {{"not/an/arg/", "error"},
       {"beta", "326"}});
    CATCH_PHRASE("is not recognized");
  }
}

