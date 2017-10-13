/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_H_
#define PAIR_H_

#include <string>
#include <vector>
#include "./space.h"
#include "./base_random.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Pair owns variables and subroutines associated with pair-wise
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
 */
class Pair : public BaseRandom {
 public:
  /// Constructor.
  /// @param rCut Interaction cut-off distance
  Pair(Space* space, const double rCut);

  /// Initialize interactions.
  //   this subroutine is written from scratch in each derived class, without
  //   any optimizations, as a test of the potential. Other potential
  //   computations and optimized and tested against this one
  virtual void initEnergy() = 0;

  /// Return total potential energy of system.
  virtual double peTot() { return peTot_; }

  // DEPRECIATED: legacy interface for old naming convention
  void Forces() { return initEnergy(); }

  /// Return potential energy of multiple particles.
  // HWH depreciate the flag?
  virtual double multiPartEner(const vector<int> multiPart, const int flag) = 0;

  /// Return potential energy of multipe particles with atom-based cutoff in 3D.
  virtual double multiPartEnerAtomCut(const vector<int> multiPart);

  /// Return potential energy of multipe particles with atom-based cutoff in 2D.
  virtual double multiPartEnerAtomCut2D(const vector<int> multiPart);

  /**
   * Compute the interaction between two particles itype and jtype separated
   * by a squared distance r2=r*r.
   * Increments the interaction in the Pair class variable "peSRone_".
   */
  virtual void multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype) {
    ASSERT(itype*jtype == r2, "multiPartEnerAtomCutInner not implemented"); }

  /// Compute the interaction between two particles itype and jtype separated
  /// by a squared distance r2=r*r. Calls multiPartEnerAtomCutInner for energy,
  /// but also computes forces using dx,dy,dz and updates fCOM_.
  virtual void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol) {
    ASSERT(0, "allPartEnerForceInner not implemented" << r2 << dx << dy
      << dz << itype << jtype << iMol << jMol); }

  /// Compute potential energy and forces of all particles.
  virtual double allPartEnerForce(const int flag) {
    ASSERT(0, "allPartEnerForce not implemented" << flag); return 0; }

  /// Compute potential energy and forces of all particles with atom-based
  /// cut-off and no cell list in 3D.
  double allPartEnerForceAtomCutNoCell();

  /// Compute potential energy and forces of all particles with atom-based
  /// cut-off and no cell list in 2D.
  double allPartEnerForceAtomCutNoCell2D();

  /// Compute potential energy and forces of all particles with atom-based
  /// cut-off and cell list in 3D.
  double allPartEnerForceAtomCutCell();

  /// Store, restore or update variables to avoid recompute
  /// of entire configuration after every trial particle move.
  /// @param uptype description of update type
  virtual void update(const vector<int> mpart,  //!< particles involved
    const int flag,   //!< type of move
    const char* uptype);

  /**
   * Update total potential energy, peTot_ with energy change, de.
   * Doesn't update anything else like neighbor lists, individual
   * contributions, etc.
   *
   * This simple pair updating method doesn't work with neighborlist. If you're
   * attempt AVB with ConfigBias, likely you forgot to enable dual-cut.
   * However, this may work with TrialGCA
   */
  virtual void update(const double de);

  /// Return total scalar virial.
  virtual double vrTot();

  /// Print XYZ format trajectory to a file.
  virtual int printxyz(const char* fileName,
    const int initFlag,  //!< open if flag is 1, append if flag is 0
    const std::string comment = "");

  /// Print GRO format trajectory to a file.
  /// @param initFlag append if 0, over-write if 1.
  virtual int printGRO(const char* fileName,
    const int initFlag = 0);

  /// Delete one particle, ipart.
  virtual void delPart(const int ipart) { delPartBase_(ipart); }

  /// Delete particles, mpart.
  virtual void delPart(const vector<int> mpart) { delPartBase_(mpart); }

  /// Add particle(s).
  virtual void addPart() { addPartBase_(); }

  /// flag to compute standard long range corrections
  // HWH: move to pair_lj.h
  int lrcFlag;

  /// Initialize pair parameters with epsilon (interaciton scale),
  /// sigma (particle size) and sigRef (reference particle size).
  void initPairParam(const vector<double> eps,
    const vector<double> sig, const vector<double> sigref);

  /// Initialize pair parameters with epsilon (interaciton scale),
  /// and sigma (particle size).
  void initPairParam(const vector<double> eps, const vector<double> sig) {
    ASSERT(sigrefFlag_ == 0, "params initialized without sigref");
    vector<double> sigref;
    initPairParam(eps, sig, sigref);
  }

  /**
   * Initialize neighbor list parameters.
   * @param neighAbove upper cut off
   * @param neighBelow lower cut off
   */
  void initNeighList(const double neighAbove, const double neighBelow);

  /// Build neighbor list for all particles.
  void buildNeighList();

  /// Return 1 if re-built neighborlist matchces current neighborlist.
  int checkNeigh();

  /**
   * Store, restore or update neighbor list variables to avoid recompute of
   * entire configuration after every particle change.
   *  flag==0, energy of old configuration (no moves, ins or dels)
   *  flag==1, new configuration, same number of particles
   *  flag==2, old configuration, preparing to delete (same as 0)
   *  flag==3, just inserted particle
   *  flag==4, old configuration, but only update the neighbor list
   *  flag==5, new configuration, but only update the neighbor list
   *
   * @param neighOneOld neighbor list to update
   */
  void updateBase(
    const vector<int> mpart,    //!< particles involved in move
    const int flag,
    const char* uptype,    //!< description of update type
    vector<vector<int> > &neigh,   //!< neighbor list to update
    vector<vector<int> > &neighOne,   //!< neighbor list to update
    vector<vector<int> > &neighOneOld);

  /// sets the cheapEnergy boolean variable
  void cheapEnergy(const int flag) {
    if (flag == 1) { cheapEnergy_ = true; } else { cheapEnergy_ = false; }; }

  /**
   * Return 1 if the currently stored energy of the configuration matches.
   *  flag=0, check currently stroed peTot_ vs recomputed with initEnergy()
   *  flag=1, check peTot_ from initEnergy() vs peTot_ from multiPartEner with
   *    mpart=all atoms
   *  flag=2, check COM forces, fCOM_ from initEnergy() vs allPartEnerForce()
   */
  int checkEnergy(const double tol, const int flag);

  /// Initialize pair data.
  void initPairData(const int natype, const vector<double> eps,
    const vector<double> sig, const vector<double> sigref);

  /// Initialize parameters with data file. Automatically searches file
  /// extensions for JSON ".json", and otherwise attempts LAMMPS data file.
  virtual void initData(const string fileName);

  /// Initialize parameters with data file. Automatically searches file
  /// extensions for JSON ".json", and otherwise attempts LAMMPS data file.
  virtual void initData(const char* fileName) { initData(string(fileName)); }

  /// Initialize parameters with LAMMPS data file.
  virtual void initLMPData(const string fileName);

  /// Initialize parameters with LAMMPS data file.
  virtual void initLMPData(const char* fileName) {
    initLMPData(string(fileName)); }

  #ifdef JSON_
    /// Initialize parameters with JSON data file.
    virtual void initJSONData(const string fileName);
    virtual void initJSONData(const char* fileName) {
      initJSONData(string(fileName)); }
  #endif  // JSON_

  /// set pair interaction parameters manually. Originally created as easy
  //  way to turn off interactions
  void epsijset(const int i, const int j, const double eps);

  /// Return excluded volume of all molecules in space, using sig_.
  double exVol(const int nGrid, const double dProbe);
  double exVol(const int nGrid) { return exVol(nGrid, 1.); }

  /// Update clusters of entire system.
  virtual void updateClusters(const double rCCut) {
    space_->peStore_ = peTot();
    space_->updateClusters(rCCut); }

  /// Initialize cut-off method by atoms or molecules.
  void initAtomCut(const int flag) {
    if (flag == 1) { atomCut_ = 1; } else { atomCut_ = 0; }; addPartBase_();
    space_->initCellAtomCut(flag);
  }

  /// Set the cut off between particles itype and jtype to rCut.
  void rCutijset(const int itype, const int jtype, const double rCut);

  /// Set type of atoms to keep track of neighbors.
  void neighTypeSet(const int itype);

  /// Initialize computation of neighborlist, neighCutOne_.
  void initNeighCut(const int flag) { neighCutOn_ = flag; }

  /// Initialize computation of potential energy map.
  void initPEMap(const int flag) { peMapOn_ = flag; }

  /**
   * Return list of molecules in neigh and their total potential energy
   * interactions. Note that peMap contains the cumulative peSRone_
   */
  void neighCutMolPEMap(vector<int> &neigh, vector<double> &peMap);

  /// Set the value of the order parameter.
  virtual void setOrder(const double order);

  /// Initialize the order parameter.
  void initOrder(
    const double order,   //!< value of the order parameter
    const char* name,     //!< name of the order parameter
    const double orderMin,  //!< minimum value of order parameter
    const double orderMax   //!< maximum value of order parameter
    ) {
    orderMin_ = orderMin; orderMax_ = orderMax; orderName_.assign(name);
    setOrder(order);
  }

  /// Read the order parameter from the comment (2nd line) of an XYZ file.
  void readXYZOrder(std::ifstream& file);

  /// Return a list of neighboring molecules of iMol, based on pair rCut.
  vector<int> iMol2neigh(const int iMol);

  /// Return the average number of neighbors within cutoff using peMap.
  double avNumNeighCut();

  /// Return a list of minimum angles between neighbors with iMol as vertex,
  /// based on pair rCut.
  vector<double> iMol2neighAngles(const int iMol);

  /** Initialize intramolecular interactions.
   *  intra_ == 0, intermolecular interactions only
   *  intra_ == 1, intramolecular interactions in addition to inter
   *  intra_ == 2, only intramolecular interactions (use with PairHybrid) */
  void initIntra(const int flag) { intra_ = flag; }
  void initIntra(const int flag, vector<vector<int> > map);
  void initIntra(vector<vector<int> > map) { initIntra(1, map); }

  /// Initialize intramolecular interactions by ignoring bonded particles.
  void initIntraBonded(const int flag);

  /// Read particle positions from XYZ file format.
  virtual void readxyz(std::ifstream& file) { space_->readxyz2(file); }

  /// Print a table file based on the potential.
  /// @param sigFac sigFac*sig is minimum separation distance
  int printTable(const char* tableFile,
    const int nElements,  //!< number of table elements
    const double sigFac = 0.9);

  /// Set all rCutij to sigmaij, which is useful for hard spheres.
  void sig2rCut();

  /// Initialize squishy tolerance when reading coordinates from a file.
  virtual void initSquishy(const int flag) { if (flag == 0) {} }

  /// Identify a particle as non physical or non physical.
  virtual void ipartNotPhysical(const int ipart) { nonphys_[ipart] = 1; }

  /// Identify a particle as non physical or non physical.
  virtual void ipartIsPhysical(const int ipart) { nonphys_[ipart] = 0; }

  /// Set all particles as physical.
  virtual void allPartPhysical() {
    std::fill(nonphys_.begin(), nonphys_.end(), 0);
  }

  /// Set flag to trigger use of a reference sigma.
  void setSigRefFlag(const int flag = 0) { sigrefFlag_ = flag; }

  // Return the pressure of the system.
  // HWH: depreciate, and may not always be up to date.
  // @param beta inverse temperature
  double pressure(const double beta);

  // read-only access to protected variables
  double f(int iAtom, int dim) const { return f_[iAtom][dim]; }  //!< force
  double rCut() const { return rCut_; }  //!< interaction cut-off distance
  double rCutMaxAll() const { return rCutMaxAll_; }
  vector<double> pe() const { return pe_; }  //!< potential energy
  vector<vector<double> > f() const { return f_; }  //!< force
  vector<vector<double> > fCOM() const { return fCOM_; }
  /// virial tensor
  vector<vector<vector<double> > > vr() const { return vr_; }
  vector<double> eps() const { return eps_; }
  vector<double> sig() const { return sig_; }
  double eps(const int i) const { return eps_[i]; }
  double sig(const int i) const { return sig_[i]; }
  double sigRef(const int i) const { return sigRef_[i]; }
  vector<vector<double> > epsij() const { return epsij_; }
  vector<vector<double> > sigij() const { return sigij_; }
  vector<vector<double> > rCutij() const { return rCutij_; }
  double rCutij(const int i, const int j) const { return rCutij_[i][j]; }
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
  Space* space() const { return space_; }

  /// Write restart file.
  virtual void writeRestart(const char* fileName) {
    writeRestartBase(fileName);
  }

  /// Constructor using restart file
  Pair(Space* space, const char* fileName);
  virtual ~Pair() {}
  virtual Pair* clone(Space* space) const = 0;

  /// reset space pointer
  virtual void reconstruct(Space* space);

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
  /// interaction spatial parameter, sigmaRef (default: unit)
  vector<double> sigRef_;
  /// interaction energy parameter, epsilon (default: unity)
  vector<vector<double> > epsij_;
  /// interaction spatial parameter, sigma (default: unity)
  vector<vector<double> > sigij_;
  /// interaction spatial parameter, sigmaRef (default: unity)
  vector<vector<double> > sigRefij_;
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

  /// Erase mpart, from neighlist, neighPtr.
  /// If atomCut_ is 1, then mpart is a list of atoms.
  /// If atomCut_ is 0, then mpart is a list of particles.
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
  int sigrefFlag_;      //!< flag to know if sigref should be read/used

  vector<int> nonphys_;  // identifies particles as non-physical, pair ignores

  /// Return true if there is an intramolecular interaction between ipart,
  /// and jpart which belong to iMol and jMol, respectively.
  bool intraCheck_(const int ipart, const int jpart,
                   const int iMol, const int jMol);

  /// defaults in constructor
  void defaultConstruction_();

  /// Delete one particle, ipart.
  void delPartBase_(const int ipart);

  /// Delete particles, mpart.
  void delPartBase_(const vector<int> mpart);

  /// Add particle(s).
  void addPartBase_();

  /// initialize neigh, neighCut and peMap.
  void initNeighCutPEMap(const vector<int> mpart);

  // Store neighCut and peMap
  // @param ii index of site in molecule
  void storeNeighCutPEMap(const int jpart,  //!< neighbor particle index
    const int ii);

  /// Write restart file.
  void writeRestartBase(const char* fileName);
};

/// Factory method.
Pair* makePair(Space* space, const char* fileName);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_H_
