/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./trial_md.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialMD::TrialMD() : Trial() {
  defaultConstruction_();
}

TrialMD::TrialMD(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
}

TrialMD::TrialMD(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
  timestep = fstod("timestep", fileName);
  nFreqZeroMomentum = fstoi("nFreqZeroMomentum", fileName);
  rescaleTemp = fstoi("rescaleTemp", fileName);
  integrator_ = fstoi("integrator", fileName);
}

void TrialMD::defaultConstruction_() {
  className_.assign("TrialMD");
  trialType_.assign("move");
  verbose_ = 0;
  timestep = 0.005;
  nFreqZeroMomentum = 1e6;
  rescaleTemp = 1;
  integrator_.assign("VelocityVerlet");
}

void TrialMD::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# timestep " << timestep << endl;
  file << "# nFreqZeroMomentum " << nFreqZeroMomentum << endl;
  file << "# rescaleTemp " << rescaleTemp << endl;
  file << "# integrator " << integrator_ << endl;
}

void TrialMD::attempt1_() {
  if (space()->nMol() <= 0) {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
    return void();
  }

  // initialize momentum, and check velocity and mass sizes
  ASSERT(vel_.size()/space()->dimen() == mass_.size(), "mass/velocity mismatch");
  if (static_cast<int>(vel_.size()) != space()->dimen()*space()->nMol()) {
    initMomentum();
  }

  // assume that forces have already been computed

  // periodically zero momentum
  if (nFreqZeroMomentum > 0) {
    if (attempted_ % nFreqZeroMomentum == 0) zeroTotalMomentum();
  }

  // integrate
  if (integrator_.compare("VelocityVerlet") == 0) {
    integrateVelocityVerlet();
  } else {
    ASSERT(0, "unrecognized integrator");
  }

  // rescale temperature
  if (rescaleTemp == 1) rescaleVelocity();

  trialAccept_();
}

void TrialMD::zeroTotalMomentum() {
  for (int dim = 0; dim < space()->dimen(); ++dim) {
    double momentum = 0;
    for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
      momentum += mass_[iMol]*vel_[space()->dimen()*iMol+dim];
    }
    for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
      vel_[space()->dimen()*iMol+dim] -= momentum
        /static_cast<double>(space()->nMol())/mass_[iMol];
    }
  }
}

double TrialMD::kineticEnergy() {
  double ke = 0.;
  for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
    for (int dim = 0; dim < space()->dimen(); ++dim) {
      const double vel = vel_[space()->dimen()*iMol+dim];
      ke += 0.5*mass_[iMol]*vel*vel;
    }
  }
  return ke;
}

double TrialMD::temperature() {
  const double dof = space()->dimen()*(space()->nMol()-1);
  return 2.*kineticEnergy()/dof;
}

void TrialMD::rescaleVelocity() {
  const double tempPrev = temperature();
  const double tempTarg = 1./criteria_->beta();
  const double factor = sqrt(tempTarg/tempPrev);
  for (int index = 0; index < static_cast<int>(vel_.size()); ++index) {
    vel_[index] *= factor;
  }
}

void TrialMD::initMomentum() {
  vel_.resize(space()->dimen()*space()->nMol(), 0.);
  mass_.resize(space()->nMol(), 1.);

  // assign velocity as uniform random number [-0.5, 0.5]
  // Note: this does not generate a gaussian distribution
  for (int index = 0; index < static_cast<int>(vel_.size()); ++index) {
    vel_[index] = uniformRanNum() - 0.5;
    // cout << "vel " << vel_[index] << endl;
  }

  zeroTotalMomentum();
  rescaleVelocity();
}

/*
 * return string for status of trial
 */
string TrialMD::printStat(const bool header) {
  stringstream stat;
  //stat << Trial::printStat(header);
  if (header) {
    stat << "ke utot temp ";
  } else {
    stat << kineticEnergy() << " "
      << kineticEnergy() + pair_->peTot() << " "
      << temperature() <<  " ";
  }
  return stat.str();
}

void TrialMD::updateFCOM() {
  if (static_cast<int>(fCOM_.size()) != space()->nMol()*space()->dimen()) {
    fCOM_.resize(space()->nMol()*space()->dimen(), 0);
  }
  std::fill(fCOM_.begin(), fCOM_.end(), 0.);
  for (int ipart = 0; ipart < space()->natom(); ++ipart) {
    const int iMol = space()->mol()[ipart];
    for (int dim = 0; dim < space()->dimen(); ++dim) {
      fCOM_[space()->dimen()*iMol + dim] += pair_->f(ipart, dim);
    }
  }
}

void TrialMD::updateVelocityHalfStep() {
  for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
    for (int dim = 0; dim < space()->dimen(); ++dim) {
      vel_[space()->dimen()*iMol + dim]
        += 0.5*timestep*fCOM(iMol, dim)*mass_[iMol];
    }
  }
}

void TrialMD::integrateVelocityVerlet() {
  // update center of mass forces
  updateFCOM();

  // obtain center of mass torques

  // update positions
  for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
    vector<double> dr(space()->dimen());
    for (int dim = 0; dim < space()->dimen(); ++dim) {
      dr[dim] = timestep*vel_[space()->dimen()*iMol + dim]
        + 0.5*timestep*timestep*mass_[iMol]*fCOM(iMol, dim);
    }
    space()->transMol(iMol, dr);
  }
  updateVelocityHalfStep();   // update half-step velocities
  pair_->initEnergy();        // compute the new forces
  updateFCOM();
  updateVelocityHalfStep();   // update half-step velocities
}

shared_ptr<TrialMD> makeTrialMD(Pair *pair, Criteria *criteria) {
  return make_shared<TrialMD>(pair, criteria);
}

shared_ptr<TrialMD> makeTrialMD() {
  return make_shared<TrialMD>();
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

