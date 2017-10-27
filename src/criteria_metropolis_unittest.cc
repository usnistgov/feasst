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
