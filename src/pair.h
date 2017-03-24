/**
 * \file
 *
 * \brief pairwise interactions
 *
 * Base class which owns variables and subroutines associated with pair-wise
 * interactions between particles.
 * The pair class calculates the pair-wise interactions between atoms using
 * an empirical, classical potential function for atoms.
 * This class owns variables and functions associated with these interactions
 * such as the forces, potential energy and virial (pressure).
 * Dimensionless units are employed such that the characteristic energy-scale
 * and length-scale of the model are unity, as well as unit mass.
 *
 * The virtual function "initEnergy()" does the main work of calculating the
 * pair-wise interactions, and is defined in the derived classes for various
 * pair-wise interaction expressions.
 * In this project, the focus is on the Lennard-Jones model, which may be
 * used to model the interaction between Argon atoms.
 */

#ifndef PAIR_H_
#define PAIR_H_

#include "./space.h"
#include "./base_all.h"

class Pair : public BaseAll {
 public:
  Pair(Space* space, const double rCut);
  Pair(Space* space, const char* fileName);
  virtual ~Pair() {}
  virtual Pair* clone(Space* space) const = 0;

  /// reset space pointer
  virtual void reconstruct(Space* space);

  /// write restart file
  virtual void writeRestart(const char* fileName) {
    writeRestartBase(fileName);
  }
  void writeRestartBase(const char* fileName);

  /// defaults in constructor
  void defaultConstruction();

  /// factory method
  Pair* makePair(Space* space, const char* fileName);

  /// pair-wise force calculation
  //   this subroutine is written from scratch in each derived class, without
  //   any optimizations, as a test of the potential. Other potential
  //   computations and optimized and tested against this one
  virtual int initEnergy() = 0;

  // DEPRECIATED: legacy interface for old naming convention
  int Forces() { return initEnergy(); }

  /// potential energy of multiple particles
  virtual double multiPartEner(const vector<int> multiPart, const int flag) = 0;
  virtual double multiPartEnerAtomCut(const vector<int> multiPart);
  virtual double multiPartEnerAtomCut2D(const vector<int> multiPart);
  virtual void multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype) {
    ASSERT(itype*jtype == r2, "multiPartEnerAtomCutInner not implemented"); }

  /// potential energy and forces of all particles
  virtual double allPartEnerForce(const int flag) {
    ASSERT(0, "allPartEnerForce not implemented" << flag); return 0; }

  /// potential energy and forces of all particles
  double allPartEnerForceAtomCutNoCell();

  /// potential energy and forces of all particles
  double allPartEnerForceAtomCutNoCell2D();

  /// potential energy and forces of all particles
  double allPartEnerForceAtomCutCell();

  /// inner loop for potential energy and forces of all particles
  virtual void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol) {
    ASSERT(0, "allPartEnerForceInner not implemented" << r2 << dx << dy
      << dz << itype << jtype << iMol << jMol); }

  /// stores, restores or updates variables to avoid order recompute
  //   of entire configuration after every change
  virtual void update(const vector<int> mpart, const int flag,
                      const char* uptype);
  virtual void update(const double de);

  /// total potential energy of system
  virtual double peTot() { return peTot_; }
  virtual double vrTot();         //!< total virial energy
  virtual int printxyz(const char* fileName, const int initFlag);
  virtual int printGRO(const char* fileName, const int initFlag);

  /// delete one particle
  void delPartBase(const int ipart);
  virtual void delPart(const int ipart) {delPartBase(ipart); }

  /// delete particles
  void delPartBase(const vector<int> mpart);
  virtual void delPart(const vector<int> mpart) {delPartBase(mpart); }

  /// add one particle
  void addPartBase();
  virtual void addPart() { addPartBase(); }

  /// flag to compute standard long range corrections
  int lrcFlag;

  /// initialize pair parameters;
  virtual void initPairParam(const vector<double> eps,
                             const vector<double> sig);

  /// full neighbor list of each molecule
  void initNeighList(const double neighAbove, const double neighBelow);
  void buildNeighList();

  int checkNeigh();  //!< check neigh list by re-building

  /// store, restore or update pair variables to avoid complete recompute
  //  after every change to the system
  void updateBase(const vector<int> mpart, const int flag, const char* uptype,
                  vector<vector<int> > &neigh, vector<vector<int> > &neighOne,
                  vector<vector<int> > &neighOneOld);

  /// sets the cheapEnergy boolean variable
  void cheapEnergy(const int flag) {
    if (flag == 1) { cheapEnergy_ = true; } else { cheapEnergy_ = false; }; }

  /// check that energy of configuration and running energy of individual
  //  particles match
  int checkEnergy(const double tol, const int flag);

  /// initialize pair data
  void initPairData(const int natype, const vector<double> eps,
                    const vector<double> sig);

  /// initialize with data file
  virtual void initData(const string fileName);
  virtual void initData(const char* fileName) { initData(string(fileName)); }

  /// initialize with LAMMPS data file
  virtual void initLMPData(const string fileName);
  virtual void initLMPData(const char* fileName) {
    initLMPData(string(fileName)); }

  /// initialize with JSON data file
  #ifdef JSON_
    virtual void initJSONData(const string fileName);
    virtual void initJSONData(const char* fileName) {
      initJSONData(string(fileName)); }
  #endif  // JSON_

  /// set pair interaction parameters manually. Originally created as easy
  //  way to turn off interactions
  void epsijset(const int i, const int j, const double eps);

  /// excluded volume of all molecules in space, using sig_
  double exVol(const int nGrid, const double dProbe);
  double exVol(const int nGrid) { return exVol(nGrid, 1.); }

  /// update clusters of entire system
  virtual void updateClusters(const double rCCut) {
    space_->updateClusters(rCCut); }

  /// initialize cut-off method by atoms or molecules
  void initAtomCut(const int flag) {
    if (flag == 1) { atomCut_ = 1; } else { atomCut_ = 0; }; addPartBase(); }

  /// set i-j cutoff
  void rCutijset(const int itype, const int jtype, const double rCut);

  /// set type of atoms to keep track of neighbors
  void neighTypeSet(const int itype);

  /// initialize computation of neighCutOne_
  void initNeighCut(const int flag) { neighCutOn_ = flag; }

  /// initialize computation of neighCutOne_
  void initPEMap(const int flag) { peMapOn_ = flag; }

  /// return list of molecules in neighCutOne and their total potential
  //  energy interactions
  void neighCutMolPEMap(vector<int> &neigh, vector<double> &peMap);

  /// set the order parameter
  virtual void setOrder(const double order);
  void initOrder(const double order, const char* name, const double orderMin,
                 const double orderMax) {
    orderMin_ = orderMin; orderMax_ = orderMax; orderName_.assign(name);
    setOrder(order);
  }

  /// initialize neigh, neighCut and peMap
  void initNeighCutPEMap(const vector<int> mpart);

  /// store neighCut and peMap
  void storeNeighCutPEMap(const int jpart, const int ii);

  /// read the order parameter from an xyz file
  void readXYZOrder(std::ifstream& file);

  /// return a vector list of neighboring molecules of iMol, based on pair rCut
  vector<int> iMol2neigh(const int iMol);

  /// compute the average number of neighbors within cutoff using peMap
  double avNumNeighCut();

  /// return a list of minimum angles between neighbors with iMol as vertex,
  //  based on pair rCut
  vector<double> iMol2neighAngles(const int iMol);

  /// initialize intramolecular interactions
  //  intra_ == 0, intermolecular interactions only
  //  intra_ == 1, intramolecular interactions in addition to inter
  //  intra_ == 2, only intramolecular interactions (use with PairHybrid)
  void initIntra(const int flag) { intra_ = flag; }
  void initIntra(const int flag, vector<vector<int> > map);
  void initIntra(vector<vector<int> > map) { initIntra(1, map); }

  /// compute intramolecular interactions (depreciated)
  double peIntra(const int iAtom);

  /// read particle positions from XYZ file format
  virtual void readxyz(std::ifstream& file) { space_->readxyz2(file); }

  /// print a table file based on the potential
  int printTable(const char* tableFile, const int nElements,
                 const double sigFac = 0.9);

  /// set all rCutij to sigmaij
  void sig2rCut();

  /// squishy tolerance when reading coordinates from a file
  virtual void initSquishy(const int flag) { if (flag == 0) {} }

  /// identify a particle as non physical or non physical
  virtual void ipartNotPhysical(const int ipart) { nonphys_[ipart] = 1; }
  virtual void ipartIsPhysical(const int ipart) { nonphys_[ipart] = 0; }
  virtual void allPartPhysical() { 
    std::fill(nonphys_.begin(), nonphys_.end(), 0);
  }

  // read-only access to protected variables
  double f(int iAtom, int dim) const { return f_[iAtom][dim]; }  //!< force
  double rCut() const { return rCut_; }  //!< interaction cut-off distance
  double rCutMaxAll() const { return rCutMaxAll_; }
  double pressure(const double beta);  //!< total pressure
  vector<double> pe() const { return pe_; }  //!< potential energy
  vector<vector<double> > f() const { return f_; }  //!< force
  vector<vector<double> > fCOM() const { return fCOM_; }
  /// virial tensor
  vector<vector<vector<double> > > vr() const { return vr_; }
  vector<double> eps() const { return eps_; }
  vector<double> sig() const { return sig_; }
  double eps(const int i) const { return eps_[i]; }
  double sig(const int i) const { return sig_[i]; }
  vector<vector<double> > epsij() const { return epsij_; }
  vector<vector<double> > sigij() const { return sigij_; }
  vector<vector<double> > rCutij() const { return rCutij_; }
  double sigij(const int i, const int j) const { return sigij_[i][j]; }
  double epsij(const int i, const int j) const { return epsij_[i][j]; }
  vector<vector<int> > neigh() const { return neigh_; }
  vector<vector<int> > neighCut() const { return neighCut_; }
  bool fastDel() const { return fastDel_; }
  int fastDelMol() const { return fastDelMol_; }
  bool neighOn() const { return neighOn_; }
  double neighAbove() const { return neighAbove_; }
  double neighBelow() const { return neighBelow_; }
  double sigMax() const { return *std::max_element(sig_.begin(), sig_.end()); }
  int atomCut() const { return atomCut_; }
  double peSRone() const { return peSRone_; }
  double order() const { return order_; }
  int orderOn() const { return orderOn_; }
  string orderName() const { return orderName_; }
  double orderMax() const { return orderMax_; }
  double orderMin() const { return orderMin_; }
  int intra() const { return intra_; }
  vector<vector<int> > contact() const { return contact_; }
  vector<int> nonphys() const { return nonphys_; }

 protected:
  Space* space_;
  int dimen_;       //!< spatial dimensionality
  double rCut_;     //!< atom interaction cut-off distance
  double rCutSq_;   //!< atom interaction cut-off distance squared
  double peTot_;    //!< total potential energy
  double peSRone_;  //!< lennard jones potential energy from subset of particles
  /// lennard jones potential energy from subset of particles
  double peSRoneAlt_;
  vector<vector<double> > f_;     //!< atomic forces
  /// center of mass force on rigid molecule
  vector<vector<double> > fCOM_;
  /// atomic potential energy
  vector<double> pe_;
  /// atomic virial tensor
  vector<vector<vector<double> > > vr_;
  /// interaction energy parameter, epsilon (default: unity)
  vector<double> eps_;
  /// interaction spatial parameter, sigma (default: unit)
  vector<double> sig_;
  /// interaction energy parameter, epsilon (default: unity)
  vector<vector<double> > epsij_;
  /// interaction spatial parameter, sigma (default: unity)
  vector<vector<double> > sigij_;
  /// potential energy cut=off for i-j type interacitons
  vector<vector<double> > rCutij_;
  /// maximum potential energy cut=off for i type interacitons
  vector<double> rCutMax_;
  double rCutMaxAll_;                 //!< maximum rCut of all types
  vector<vector<int> > neighCut_;      //!< list of neighbors for each particle
  /// list of neighbors for subset of particles
  vector<vector<int> > neighCutOne_;
  /// list of neighbors for subset of particles
  vector<vector<int> > neighCutOneOld_;
  vector<vector<int> > neigh_;      //!< list of neighbors for each particle
  /// list of neighbors for subset of particles
  vector<vector<int> > neighOne_;
  /// list of neighbors for subset of particles
  vector<vector<int> > neighOneOld_;
  bool neighOn_;                                    //!< is neighbor list on
  int neighCutOn_;                      //!< flag to compute neighCut
  int peMapOn_;                         //!< flag to compute peMap
  vector<vector<double> > peMap_;          //!< potential energy interactions
  double neighAbove_;    //!< distance cut off distance for neighbor criteria
  /// squared distance cut off distance for neighbor criteria
  double neighAboveSq_;
  /// distance cut off distance for neighbor criteria
  double neighBelow_;
  /// squared distance cut off distance for neighbor criteria
  double neighBelowSq_;
  /// neighbor list only keeps track of ij interactions.
  //  If empty, all interactions are considered
  vector<int> neighType_;
  int neighTypeScreen_;               //!< screen neighbors by type if == 1
  
  /// erase molecule from neighlist
  void eraseNeigh_(const vector<int> mpart, vector<vector<int> > *neighPtr);
  
  bool fastDel_;                   //!< use fast method of deleting particles
  int fastDelMol_;                   //!< molecule last deleted by fast method
  /// cheap potential calculation for CBMC multiple first-bead insertions
  bool cheapEnergy_;
  /// atomic or molecule distance based cut-off, neigh list is by atom if 1,
  //  or by molecule if 0
  int atomCut_;
  /// list of files used to initialize with data files
  vector<string> initDataFiles_;
  /// record list of epsijset for restart file
  vector<vector<double> > epsijsetRecord_;
  
  int intra_;         //!< flag for intramolecular interactions

  // order parameter variables
  double order_, orderMin_, orderMax_;
  int orderOn_;     //!< toggle order parameter
  string orderName_;

  // contact map
  vector<vector<int> > contact_;    //!< contact map between pairs of particles
  /// contact pbc map between pairs of particles
  vector<vector<vector<double> > > contactpbc_;

  int sigepsDefault_;   //!< flag to know if default eps, sig are assigned

  vector<int> nonphys_;  // identifies particles as non-physical, pair ignores

  /// check if there is an intramolecular interaction that is allowed
  bool intraCheck_(const int ipart, const int jpart,
                   const int iMol, const int jMol);

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);}
};

#endif  // PAIR_H_

