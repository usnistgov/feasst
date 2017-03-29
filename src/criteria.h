/**
 * \file
 *
 * \brief acceptance criteria for monte carlo trials
 *
 */

#ifndef CRITERIA_H_
#define CRITERIA_H_

#include "./base_all.h"

class Space;
class Pair;

class Criteria : public BaseAll {
 public:
  Criteria(const double beta, const double activ);
  explicit Criteria(const char* fileName);
  virtual ~Criteria() {}
  virtual Criteria* clone() const = 0;

  /// write restart file
  virtual void writeRestart(const char* fileName) {
    writeRestartBase(fileName);
  }
  void writeRestartBase(const char* fileName);

  /// factory method
  Criteria* makeCriteria(const char* fileName);

  /// defaults in constructor
  void defaultConstruction();

  /// acceptance criteria for trial moves
  virtual int accept(const double lnpMet, const double peNew,
    const char* moveType, const int reject) = 0;

  /// store macrostate variables of old configuration
  virtual void store(const Space* space, Pair* pair);

  /// zero all statistics and accumulators
  virtual void zeroStat() {}

  void betaset(const double beta) { beta_ = beta; }
  void activset(const double activ) { activ_ = activ; }
  void pressureset(const double pressure) {
    pressureFlag_ = 1; pressure_ = pressure;
  }

  /// activity for each molecule type
  vector<double> activVec;
  void addActivity(const double activ) {
    activVec.push_back(activ);
  }

  /// read only access to protected variables
  double beta() const { return beta_; }
  double activ() const;
  double activ(const int type) const { return activVec[type]; }
  double pressure() const { return pressure_; }
  int pressureFlag() const { return pressureFlag_; }
  int printBeta() const { return printBeta_; }
  int printPressure() const { return printPressure_; }

 protected:
  double beta_;
  double activ_;      //!< exp(beta*mu)/lambda^3
  double pressure_;
  int pressureFlag_;

  // flag to print in log files
  int printBeta_;
  int printPressure_;

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);}

  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const = 0;
};

#endif  // CRITERIA_H_

