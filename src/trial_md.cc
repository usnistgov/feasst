#include "./trial_md.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialMD::TrialMD() : Trial() {
  defaultConstruction();
}
TrialMD::TrialMD(Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria) {
  defaultConstruction();
}
TrialMD::TrialMD(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  defaultConstruction();
  timestep = fstod("timestep", fileName);
  nFreqZeroMomentum = fstoi("nFreqZeroMomentum", fileName);
  rescaleTemp = fstoi("rescaleTemp", fileName);
  integrator_ = fstoi("integrator", fileName);
}

/**
 * default construction
 */
void TrialMD::defaultConstruction() {
  className_.assign("TrialMD");
  trialType_.assign("move");
  verbose_ = 0;
  timestep = 0.005;
  nFreqZeroMomentum = 1e6;
  rescaleTemp = 1;
  integrator_.assign("VelocityVerlet");
}

/**
 * write restart file
 */
void TrialMD::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# timestep " << timestep << endl;
  file << "# nFreqZeroMomentum " << nFreqZeroMomentum << endl;
  file << "# rescaleTemp " << rescaleTemp << endl;
  file << "# integrator " << integrator_ << endl;
}

/**
 * Attempt trial
 */
void TrialMD::attempt1() {
  if (space_->nMol() <= 0) {
    trialMoveDecide(0, 0);   // ensured rejection, however, criteria can update
    return void();
  }

  // initialize momentum, and check velocity and mass sizes
  ASSERT(vel_.size()/space_->dimen() == mass_.size(), "mass/velocity mismatch");
  if (static_cast<int>(vel_.size()) != space_->dimen()*space_->nMol()) {
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

  trialAccept();
}

/**
 * zero total momentum of the system
 */
void TrialMD::zeroTotalMomentum() {
  for (int dim = 0; dim < space_->dimen(); ++dim) {
    double momentum = 0;
    for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
      momentum += mass_[iMol]*vel_[space_->dimen()*iMol+dim];
    }
    for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
      vel_[space_->dimen()*iMol+dim] -= momentum
        /static_cast<double>(space_->nMol())/mass_[iMol];
    }
  }
}

/**
 * compute the kinetic energy
 */
double TrialMD::kineticEnergy() {
  double ke = 0.;
  for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
    for (int dim = 0; dim < space_->dimen(); ++dim) {
      const double vel = vel_[space_->dimen()*iMol+dim];
      ke += 0.5*mass_[iMol]*vel*vel;
    }
  }
  return ke;
}

/**
 * compute the temperature
 * NOTE: This assumes monotonic, isotropic system for number of degrees of
 * freedom.
 */
double TrialMD::temperature() {
  const double dof = space_->dimen()*(space_->nMol()-1);
  return 2.*kineticEnergy()/dof;
}

/**
 * rescale velocity according to the temperature
 */
void TrialMD::rescaleVelocity() {
  const double tempPrev = temperature();
  const double tempTarg = 1./criteria_->beta();
  const double factor = sqrt(tempTarg/tempPrev);
  for (int index = 0; index < static_cast<int>(vel_.size()); ++index) {
    vel_[index] *= factor;
  }
}

/**
 * initialize velocity, the lazy way
 * Maxwell-Boltzmann is expected to rapidly estability after a few steps
 */
void TrialMD::initMomentum() {
  vel_.resize(space_->dimen()*space_->nMol(), 0.);
  mass_.resize(space_->nMol(), 1.);

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

/**
 * update center of mass forces
 */
void TrialMD::updateFCOM() {
  if (static_cast<int>(fCOM_.size()) != space_->nMol()*space_->dimen()) {
    fCOM_.resize(space_->nMol()*space_->dimen(), 0);
  }
  std::fill(fCOM_.begin(), fCOM_.end(), 0.);
  for (int ipart = 0; ipart < space_->natom(); ++ipart) {
    const int iMol = space_->mol()[ipart];
    for (int dim = 0; dim < space_->dimen(); ++dim) {
      fCOM_[space_->dimen()*iMol + dim] += pair_->f(ipart, dim);
    }
  }
}

/**
 * update half-step velocities
 */
void TrialMD::updateVelocityHalfStep() {
  for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
    for (int dim = 0; dim < space_->dimen(); ++dim) {
      vel_[space_->dimen()*iMol + dim]
        += 0.5*timestep*fCOM(iMol, dim)*mass_[iMol];
    }
  }
}

/**
 * velocity verlet integration
 */
void TrialMD::integrateVelocityVerlet() {
  // update center of mass forces
  updateFCOM();

  // obtain center of mass torques

  // update positions
  for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
    vector<double> dr(space_->dimen());
    for (int dim = 0; dim < space_->dimen(); ++dim) {
      dr[dim] = timestep*vel_[space_->dimen()*iMol + dim]
        + 0.5*timestep*timestep*mass_[iMol]*fCOM(iMol, dim);
    }
    space_->transMol(iMol, dr);
  }
  updateVelocityHalfStep();   // update half-step velocities
  pair_->initEnergy();        // compute the new forces
  updateFCOM();
  updateVelocityHalfStep();   // update half-step velocities
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

