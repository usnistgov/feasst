/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef PAIR_TABULAR_H_
#define PAIR_TABULAR_H_

#include "./pair.h"
#include "./table.h"

namespace feasst {

class PairTabular : public Pair {
 public:
  explicit PairTabular(Space* space);
  PairTabular(Space* space, const char* fileName);
  virtual ~PairTabular() {}
  virtual PairTabular* clone(Space* space) const {
    PairTabular* p = new PairTabular(*this);
    p->reconstruct(space);
    return p;
  }

  virtual int initEnergy() { return 0; }

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag) {
    if (flag == 0 || multiPart.size() < 2) {} return 0;
  }

  /// initialize tables
  void initTableHard(const int iType, const int jType, const char* fileName);
  void initTableCut(const int iType, const int jType, const char* fileName);
  void initTablePE(const int iType, const int jType, const char* fileName);

  /// initialize whether to include more than hard interactions
  void initHard(const int iType, const int jType, const int flag = 1);
  void initHard(const int flag = 1) { initHard(0, 0, flag); }

  /// read particle positions from XYZ file format
  void readxyzeuler(std::ifstream& file) {}
  
  /// squishy tolerance when reading coordinates from a file
  void initSquishy(const int flag) { squishy_ = flag; }

  /// compute the second virial coefficient directly from the table
  double b2(const double beta, const int expandt = 1, const int expandz = 1) { return beta*expandt*expandz; }

  /// initialize the potential cutoffs
  void initCuts() {}

  /// convert particle type to table ID via tabID_
  int type2tabID(const int iType) { return iType; }

  /// read-only access of protected variables
  vector<vector<int> > hardFlag() const { return hardFlag_; }
  vector<vector<double> > rCutInner() const { return rCutInner_; }
  vector<vector<shared_ptr<Table> > > tabHard() const { return tabHard_; }
  vector<vector<double> > peMin() const { return peMin_; }
  vector<vector<string> > tabType() const { return tabType_; }

 protected:
  vector<vector<double> > rCutInner_;     //!< hard sphere below inner cut-off
  double deSR_;             //!< lennard jones potential energy change
  int squishy_;             //!< imprecise coordinates form file
  /// table which defines hard contact distance as a function of orientation
  vector<vector<shared_ptr<Table> > > tabHard_;
  /// table which defines cut-off distance as a function of orientation
  vector<vector<shared_ptr<Table> > > tabCut_;
  /// tablef or potential energy as function of z=(r-rh)/(rc-rh)
  vector<vector<shared_ptr<Table> > > tabPE_;
  vector<vector<string> > tabHardFileName_;
  vector<vector<string> > tabCutFileName_;
  vector<vector<string> > tabPEFileName_;
  vector<vector<int> > hardFlag_;              //!< no depletant
  vector<vector<string> > tabType_;        //!< type of table
  vector<vector<double> > peMin_;    //!< minimum of potential energy
  vector<int> tabID_;             //!< table index to site type
  double tol_;      //!< table tolerance
};

}  // namespace feasst

#endif  // PAIR_TABULAR_H_

