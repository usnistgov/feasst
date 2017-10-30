/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * constants from CODATA
 * http://dx.doi.org/10.1103/RevModPhys.84.1527
 * http://dx.doi.org/10.1063/1.4724320
 */
class PhysicalConstants {};  // here is where we trick ../tools/makeFactory.sh.

const double boltzmannConstant = 1.3806488E-23; //!< J/K
const double avogadroConstant = 6.02214129E+23; //!< 1/mol
const double elementaryCharge = 1.602176565E-19;  //!< C
const double permitivityVacuum = 8.854187817E-12;  //!< C^2/J/m
const double idealGasConstant = boltzmannConstant*avogadroConstant;  //!< J/K/mol
const double joulesPercal = 4.184;  //!< J/cal

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  //PHYSICAL_CONSTANTS_H_
