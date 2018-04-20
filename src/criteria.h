/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef CRITERIA_H_
#define CRITERIA_H_

#include "./functions.h"
#include "./base_random.h"
#include "./space.h"
#include "./pair.h"

namespace feasst {

/**
 * Acceptance criteria for Monte Carlo trials.
 */
class Criteria : public BaseRandom {
 public:

  // Constructor
  // HWH: Depreciate in favor of args
  Criteria(
    /// inverse temperature: \f$ \beta = \frac{1}{k_B T} \f$
    const double beta,
    /** activity: \f$ z = \frac{exp(\beta \mu)}{\Lambda^3} \f$.
     *  Only use this activity with one component fluctuating in number. */
    const double activ);

  /// Constructor
  explicit Criteria(
    /// inverse temperature: \f$ \beta = \frac{1}{k_B T} \f$
    const double beta,
    /**
     * allowed string key pairs (e.g., dictionary)
     *
     * activ (or activ0) : activity of the first particle type
     *   \f$ z = \frac{exp(\beta \mu)}{\Lambda^3} \f$.
     */
    const argtype &args = argtype());

  /// Return inverse temperature.
  double beta() const { return beta_; }

  /// Set the inverse temperature.
  void betaset(const double beta) { beta_ = beta; }

  /** Return activity.
   *  Only use this activity with one component fluctuating in number. */
  // NOTE HWH: remove activity from constructor. Its no longer a single
  // quantity. Or use args?
  double activ() const;

  /** Set the activity.
   *  Only use this activity with one component fluctuating in number. */
  void activset(const double activ) { activ_ = activ; }

  /// Add an activity in order of type.
  void addActivity(const double activ) {
    activVec_.push_back(activ);
  }

  /// Return activity of a given type, where type is index in order of add.
  double activ(const int type) const;

  /// Return the number of activies that were added
  int nActiv() const { return static_cast<int>(activVec_.size()); }

  /** Return whether to accept (1) or reject (0) the trial based randomly on
   *  the acceptance criteria. */
  virtual int accept(
    /// logarithm of the Metropolis acceptance probability
    const double lnpMet,
    const double peNew,  //!< potential energy of proposed configuration
    const char* moveType,
    const int reject  //!< automatically reject if == 1
    ) = 0;

  /// Construct by checkpoint file.
  explicit Criteria(const char* fileName);

  virtual ~Criteria() {}
  virtual Criteria* clone() const = 0;

  /// Write restart file.
  virtual void writeRestart(const char* fileName) {
    writeRestartBase(fileName);
  }
  void writeRestartBase(const char* fileName);

  /// Store macrostate variables of old configuration.
  virtual void store(Pair* pair);

  /// Zero all statistics and accumulators.
  virtual void zeroStat() {}

  /// Set the pressure.
  void pressureset(const double pressure) {
    pressureFlag_ = 1; pressure_ = pressure;
  }

  /// Return the pressure.
  double pressure() const { return pressure_; }

  /// Return 1 if pressure has been set.
  int pressureFlag() const { return pressureFlag_; }

 protected:
  double beta_;
  double activ_;
  double pressure_;
  int pressureFlag_;

  /// Activity for each molecule type.
  vector<double> activVec_;

  /// defaults in constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const = 0;
};

/// Factory method.
Criteria* makeCriteria(const char* fileName);

// parse beta from args
double parseBeta_(argtype * args);

}  // namespace feasst

#endif  // CRITERIA_H_

