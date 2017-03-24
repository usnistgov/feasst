/**
 * \file
 *
 * \brief trial add particles for Monte Carlo
 *
 */

#ifndef TRIAL_CONFIGBIAS_H_
#define TRIAL_CONFIGBIAS_H_

#include "./trial.h"

class Space;
class Pair;
class Criteria;
class MC;

class TrialConfigBias : public Trial {
 public:
  explicit TrialConfigBias(const char* molType);
  TrialConfigBias(Space *space, Pair *pair, Criteria *criteria,
                  const char* molType);
  TrialConfigBias(const char* fileName, Space *space, Pair *pair,
                  Criteria *criteria);
  ~TrialConfigBias() {}
  TrialConfigBias* clone(Space* space, Pair *pair, Criteria *criteria) const;
  shared_ptr<TrialConfigBias> cloneShrPtr
    (Space* space, Pair* pair, Criteria* criteria) const;

  /// attempt insertion
  void attempt1();

  /// implement one step in cbmc, return list of atoms that were moved
  vector<int> buildOneStep(const int step, const int flag);

  /// implement all steps in cbmc, return energy and rosenbluth factor
  void buildAllStep(double &pe, double &lnw, const int newConf);

  /// type of grow move
  string growType;

  /// add generic growth step
  void addStepGeneric(const int iAtom, const int n);

  /// add growth step, placing iAtom randomly in box n times
  void addStepBox(const int iAtom, const int n);

  /// add growth step, placing iAtom along bond with jAtom, n times
  void addStepBond(const int iAtom, const int jAtom, const int n);

  /// add growth step, placing iAtom along angle with jAtom and kAtom, ntimes
  void addStepAngle(const int iAtom, const int jAtom,
                    const int kAtom, const int n);

  /// add growth step, placing iAtom and jAtom along branch
  //  with kAtom and lAtom, ntimes
  void addStepBranch(const int iAtom, const int jAtom, const int kAtom,
                     const int lAtom, const int n);

  /// add growth step, placing iAtom at COM for free
  void addStepFreeCOM(const int iAtom);

  /// add aggregation volume bias (avb) step, placing iAtom in the avb of the
  // list of target atoms, tmpart, with avb defined by rAove and rBelow, n times
  void addStepAVB(const int iAtom, const vector<int> tmpart,
                  const double rAbove, const double rBelow, const int n);

  /// update stepParam, using current values from space class
  void updateStepParams();

  /// reduce list of particles for a molecule to only particles
  //  involved in CB move
  void mpartCBparse();

  /// initialize dual-cut configurational bias
  void initDualCut(const int flag);

  /// initialize trials in given MC class based on reading input file
  void file2MC(const char* fileName, MC* mc);
  void lmp2MC(const char* fileName, MC* mc);
  void json2MC(const char* fileName, MC* mc);

  /// flag on how to handle insert/delete trials
  //   ==0, ignore insert/delete trials
  //   ==1, add insert/delete trials
  //   ==2, convert insert/delete to regrow
  int insDelFlag;

  /// exponent for anglar or bond intramolecular interaction
  //  U=k(x-x0)^exponent
  int bondExponent;
  int angleExponent;
  
  /// read only access to protected variables
  string molType() const { return molType_; }
  int nSteps() const { return static_cast<int>(stepAtom_.size()); }
  string stepName(const int i) const { return stepName_[i]; }

 protected:
  string molType_;              //!< type of molecule to add
  vector<int> stepAtom_;        //!< atom to place in each step
  vector<int> stepTrials_;      //!< number of trials to preform for each step
  vector<string> stepName_;     //!< list of configurational bias step names
  vector<vector<double> > stepParam_;   //!< list of configurational bias
                                        // step parameters, for each step
  vector<vector<int> > stepList_;   //!< list of atoms involved in step
  vector<int> mpartall_;            //!< all particles in selected molecule
  int dualCut_;       //!< dual-cut Configurational Bias
  string regionOld_;  //!< AVB region for old configuration
  vector<int> tmpartOld_;   //!< target particle for old configuration
  double pairOrderOld_;   //!< track order parameter in pair class to see
                          // if stepParm needs to be updated

  /// for each step, list the number of consequentive steps to perform
  //   consequentive steps have are determined by 1 trial moves in data file
  vector<int> stepOrder_;

  /// initialize buildAllStep for AVB insertions/deletions, return prefactor
  double buildAVB(const int newConf);

  /// initialize buildAllStep for AVB2 moves, return prefactor
  double buildAVB2();

  /// initialize buildAllStep for AVB3 moves, return prefactor
  double buildAVB3();

  /// convert tmpart from list of all sites in a target molecule to a random
  //  site available for AVB. Also return neighbors of that same type
  vector<int> mpart2neigh(vector<int> &tmpart);

  /// randomly choose an atom in list to define mpart
  void chooseMPart(const vector<int> &list);

  /// return list of atoms in the out region mpart
  vector<int> iAtom2notNeigh(const int iAtom, const vector<int> &neigh);

  /// clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const;
};

#endif  // TRIAL_CONFIGBIAS_H_

