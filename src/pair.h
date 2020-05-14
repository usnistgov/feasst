/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_H_
#define PAIR_H_

#include <string>
#include <vector>
#include "./space.h"
#include "./base_random.h"
#include "./accumulator.h"

namespace feasst {

/**
 * Pair owns variables and subroutines associated with pair-wise
 * interactions between particles.
 * The pair class calculates the pair-wise interactions between atoms using
 * an empirical, classical potential function for atoms.
 * This class owns variables and functions associated with these interactions
 * such as the forces, potential energy and virial (pressure).
 *
 * All Pair classes are expected to have the following
 *  - initialization via data file, which sets the epsilons and sigmas
 *  - interaction cutoff distance
 *
 * In order to create a custom Pair, you can follow two similar procedures.
 *
 * First, you may define a custom Pair code in the same file as
 * "int main()" for C++. See an example of this for the Jagla potential
 * <a href="tutorial/2_jagla_1_ref-config_README.html">test case</a>.
 *
 * Second, you may copy existing "pair_*" files and replace the class name
 * and header guards (e.g. BASECLASS_DERIVED_H_).
 */
class Pair : public BaseRandom {
 public:
  /// Constructor.
  Pair(Space* space,
    /**
     * allowed string key pairs (e.g. dictionary):
     *
     * rCut : floating point value of the interaction truncation distance.
     *        Typically used by only the simplest pure-component models.
     *
     *  - (default): 0.
     */
    const argtype &args = argtype());

  /// Initialize interactions.
  //   this subroutine is written from scratch in each derived class, without
  //   any optimizations, as a test of the potential. Other potential
  //   computations and optimized and tested against this one
  virtual void initEnergy() { peTot_ = allPartEnerForce(2); }

  /// Return total potential energy of system.
  virtual double peTot() { return peTot_; }

  // DEPRECIATED: legacy interface for old naming convention
  void Forces() { return initEnergy(); }

  /// Return potential energy of multiple particles.
  virtual double multiPartEner(const vector<int> multiPart, const int flag);

  /**
   * Compute the interaction between two particles itype and jtype separated
   * by a squared distance r2=r*r.
   * Increments the interaction in the Pair class variable "peSRone_".
   * Return 1 if particles are neighbors.
   */
  virtual int multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype) {
    ASSERT(itype*jtype == r2, "multiPartEnerAtomCutInner not implemented");
    return 1; }

  /// Compute the interaction between two particles itype and jtype separated
  /// by a squared distance r2=r*r. Calls multiPartEnerAtomCutInner for energy,
  /// but also computes forces using dx,dy,dz and updates fCOM_.
  virtual void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol) {
    if (dx*dy*dz*iMol*jMol == 0) {}  // remove unused variable warning
    multiPartEnerAtomCutInner(r2, itype, jtype); }

  /// Compute potential energy and forces of all particles.
  virtual double allPartEnerForce(
    /// If flag == 0, no compute. return the stored potential energy (peTot_)
    /// If flag == 1, compute
    /// If flag == 2, compute without optimizations (e.g., cells, etc)
    const int flag = 1);

  /// Compute potential energy and forces of all particles with molecule-based
  /// cut-off and no cell list in 3D.
  double allPartEnerForceNoCell();

  /// Compute interactions between two molecules.
  /// For the default behavior in Pair, loop through each site and consider
  /// rCut between each site.
  virtual void allPartEnerForceMolCutInner(const double r2,
    const int iMol, const int jMol, const double dx,
    const double dy, const double dz);

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

  // For CriteriaMayer, simply set the total potential energy.
  virtual void updatePeTot(const double peTot);

  // Return total scalar virial.
  // HWH: Depreciate or refactor.
  virtual double vrTot();

  /// Print XYZ format trajectory to a file.
  virtual int printXYZ(const char* fileName,
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
  virtual void addMol(const char* fileName) {
    space_->addMol(fileName); addPartBase_(); }
  virtual void addPart() { addPartBase_(); }
  void addMol(const std::string fileName) { addMol(fileName.c_str()); }
  virtual void addMol() { space_->addMol(); addPartBase_(); }

  /// Add molecule at the specified position. If no file provided, use Space
  /// default behavior to decide which molecule (typically the first).
  void addMol(const vector<double> position, const char* fileName = "");

  /// Add molecule at the specified position, using std::string.
  void addMol(const vector<double> position, const std::string fileName) {
    addMol(position, fileName.c_str()); }

  /// flag to compute standard long range corrections
  // HWH: move to pair_lj.h, but this requires that PairLJCoulEwald inherits it
  int lrcFlag;

  /// Initialize pair parameters with epsilon (interaciton scale),
  /// sigma (particle size) and sigRef (reference particle size).
  virtual void initPairParam(const vector<double> eps,
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
  virtual void cheapEnergy(const int flag) {
    if (flag == 1) { cheapEnergy_ = true; } else { cheapEnergy_ = false; }; }

  /**
   * Return 1 if the currently stored energy of the configuration matches.
   *  flag=0, check currently stroed peTot_ vs recomputed with initEnergy()
   *  flag=1, check peTot_ from initEnergy() vs peTot_ from multiPartEner with
   *    mpart=all atoms. This fails for intramolecular interactions.
   *  flag=2, check COM forces, fCOM_ from initEnergy() vs allPartEnerForce()
   */
  int checkEnergy(const double tol, const int flag);

  /// Initialize pair data.
  void initPairData(const int natype, const vector<double> eps,
    const vector<double> sig, const vector<double> sigref);

  /// Initialize parameters with data file. Automatically searches file
  /// extensions for JSON ".json", and otherwise attempts LAMMPS data file.
  /// Also uses the data file to initialize the Space class.
  virtual void initData(const string fileName);

  /// Initialize parameters with data file. Automatically searches file
  /// extensions for JSON ".json", and otherwise attempts LAMMPS data file.
  /// Also uses the data file to initialize the Space class.
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

  /// Set pair interaction parameters manually. Originally created as easy
  /// way to turn off interactions.
  /// Note: The site types in LMP start from 1, but FEASST starts from 0.
  void epsijset(const int iSiteType, const int jSiteType, const double eps);

  // Set the pair interaction parameters manually. Originally created for field
  // Note: may have trouble with restarts.
  void initEps(const int siteType, const double eps);
  void initSig(const int siteType, const double sig);

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
    space_->initAtomCut(flag);
  }

  /// Set the cut off between particles itype and jtype to rCut.
  void rCutijset(const int itype, const int jtype, const double rCut);

  /// Set the cut off between particles of itype to rCut.
  void initRCut(const int itype, const double rCut) {
    rCutijset(itype, itype, rCut); }

  /// Set the cut off between all particles types equal.
  void equateRcutForAllTypes();

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
  virtual void readXYZ(std::ifstream& file) {
    space_->readXYZ(file);
    addPart(); }

  /// Scale domain (see Space)
  virtual void scaleDomain(const double factor, const int dim) {
    space()->scaleDomain(factor, dim); }

  /// Scale domain (see Space)
  virtual void scaleDomain(const double factor) { space()->scaleDomain(factor); }

  /// Print a table file based on the potential.
  /// @param sigFac sigFac*sig is minimum separation distance
  int printTable(const char* tableFile,
    const int nElements,  //!< number of table elements
    const double sigFac = 0.9);

  /// Set all rCutij to sigmaij, which is useful for hard spheres.
  void sig2rCut();
  void halfSig2rCut();

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

  // Accumulate averages
  virtual void accumulate();
  virtual void accumulateReset();
  Accumulator accumulator() const { return peAcc_; }

  // read-only access to protected variables
  double f(int iAtom, int dim) const { return f_[iAtom][dim]; }  //!< force
  double rCut() const { return rCut_; }  //!< interaction cut-off distance
  double rCutMaxAll() const { return rCutMaxAll_; }
  vector<double> pe() const { return pe_; }  //!< potential energy
  vector<vector<double> > f() const { return f_; }  //!< force
  vector<vector<double> > fCOM() const { return fCOM_; }
  // virial tensor: HWH: depreciate or refactor
  vector<vector<vector<double> > > vr() const { return vr_; }
  vector<double> eps() const { return eps_; }
  vector<double> sig() const { return sig_; }
  double eps(const int i) const { return eps_[i]; }
  double sig(const int i) const { return sig_[i]; }
  double sigRef(const int i) const { return sigRef_[i]; }
  vector<vector<double> > epsij() const { return epsij_; }
  vector<vector<double> > sigij() const { return sigij_; }
  vector<vector<double> > sigRefij() const { return sigRefij_; }
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
  virtual Pair* clone(Space* space) const {
    if (space == nullptr) {};  // remove unused parameter warning
    ASSERT(0, "clone not implemented");
    return nullptr; }

  /// reset space pointer
  virtual void reconstruct(Space* space);

  /// Compute forces if 1 (default 0)
  virtual void initForces(const int flag) { forcesFlag_ = flag; }

  // potential energy interactions
  vector<vector<double> > peMap() const { return peMap_; }
  vector<vector<int> > neighCutOne() const { return neighCutOne_; }

  // if flag == 0, dont skip interaction if epsilon == 0
  // if flag == 1 (default), skip if epsilon == 0
  void setSkipEPS0(const int flag) { skipEPS0_ = flag; }

 protected:
  Space* space_;
  int dimen_;       //!< spatial dimensionality
  double rCut_;     //!< atom interaction cut-off distance
  double rCutSq_;   //!< atom interaction cut-off distance squared
  double peTot_;    //!< total potential energy
  double peSR_;
  double deSR_;
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
  /// record list of initEps for restart file
  vector<vector<double> > epsRecord_;
  vector<vector<double> > sigRecord_;

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

  /// Flag which tells pair not to skip count eps=0
  int skipEPS0_ = 1;

  /// VMD labels read for data file
  vector<string> VMDlabels_;

  /*
   Here begins the pair-wise "loop" routines for interactions between particles
   and sites. Note, particles may also be referred to as molecules.

   The different types of loops are listed below, with their nomenclature:
   1. Site vs Mol :: site-based vs molecule-based loop

   If a site has eps == 0 and skipEPS0 != 0, then it is skipped.
   Attempts to use cell list and/or neighborlist if available
   */
  virtual double pairLoopSite_(
    /// must be sorted (due to set_difference operation)
    const vector<int> &siteList,
    /// set to 1 to force no use of cell list (error checking)
    const int noCell = 0);

  /// Overload default such that siteList is all sites in space
  double pairLoopSite_(const int noCell = 0) {
    const vector<int> &allSites = space_->listAtoms();
    return pairLoopSite_(allSites, noCell); }

  /**
   * See pairLoopSite_, except the cut-off is based on particles (e.g., molecules).
   * Particle cut-off is measured from the first site in a particle.
   * The siteList is assumed to be continguous in molecule number.
   * HWH NOTE: peMap is not implemented.
   */
  double pairLoopParticle_(const vector<int> &siteList, const int noCell = 0);

  /// Compute the interaction between two sites
  virtual void pairSiteSite_(
    const int &iSiteType,  //!< type of first site
    const int &jSiteType,  //!< type of second site
    double * energy,      //!< energy of interaction
    double * force,       //!< force/rij of interaction
    int * neighbor,       //!< 1 if neighbor, 0 otherwise
    const double &dx,      //!< x-dimension separation
    const double &dy,      //!< y-dimension separation
    const double &dz       //!< z-dimension separation
    ) { ASSERT(0, "not implemented"); }

  /// Return whether or not to use a cell list
  bool useCellForSite_(
    /// Particle type. Use -1 for all types.
    const int itype = -1);

  /// Tag a site as a neighbor of another site within Subset loop
  void setNeighbor_(const double &r2,   //!< squared separation distance
    const int &siteIndex,  //!< index of siteList in pairSubset
    const int &neighSite,   //!< (space) site index of neighboring site
    const int &siteType,    //!< type of site in siteList
    const int &neighType    //!< type of site for neighbor
    );

  int forcesFlag_ = 0;  // compute forces if == 1

  // accumulator for average energy from MC trials
  Accumulator peAcc_;

  // Cheaply compute the interaction between two particles (e.g., dual cut)
  virtual void pairParticleParticleCheapEnergy_(const double &r2, const int &itype,
    const int &jtype, double * energy, double * force) {
    ASSERT(itype*jtype*0 + 1 == (r2 + *energy + *force)*0 + 2,
           "not implemented"); }
};

/// Factory method with implementation generated at build time.
Pair* makePair(Space* space, const char* fileName);

}  // namespace feasst

#endif  // PAIR_H_
