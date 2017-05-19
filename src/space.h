#ifndef SRC_SPACE_H_
#define SRC_SPACE_H_

#include <vector>
#include <memory>
#include <string>
#include "./functions.h"
#include "./histogram.h"
#include "./base_all.h"
#ifdef XDRFILE_H_
  extern "C" {
    #include "xdrfile.h"
    #include "xdrfile_xtc.h"
    #include "xdrfile_trr.h"
  }
#endif  // XDRFILE_H_

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * The space class owns variables and functions associated with the real-space
 * position of particles and the domain in which they reside.
 */
class Space : public BaseAll {
 public:
  /// Construct with spatial dimension and ID number.
  Space(int dimen = 3, int id = 0);

  /// Construct by checkpoint file.
  explicit Space(const char* fileName);
  ~Space();

  /** Return a deep copy of self. Note that this pointer was constructed with
   *  the "new" directive, which means that it requires a subsequent delete to
   *  avoid a memory leak. */
  Space* clone() const;

  /** Return a deep copy of self as shared pointer with automated memory
   *  management (e.g., preferred over the clone method). */
  shared_ptr<Space> cloneShrPtr() const;

  /** Write restart file. Print each molecule type that was added
   *  followed by the number of that kind of molecule.
   *  Finally, print coordinates of all molecules in that order. */
  void writeRestart(const char* fileName);

  /* NOTE: HWH: Depreciated function. Don't use this.
   * Initialize the atomic positions according to some formula used for
   * testing purposes. Use python for more general input. */
  int init_config(const int natom);

  /// Write the atomic positions to a file.
  int printXYZ(const char* fileName,
    const int initFlag,
   /*  Open for first time and write vmd script if flag is 1.
    *  Append if flag is 0.
    *  Print vmd file for xtc if flag is 2.
    */
    const std::string comment="");

  // NOTE HWH: Depreciated for style. Use printXYZ.
  int printxyz(const char* fileName, const int initFlag,
    const std::string comment="") { return printXYZ(fileName, initFlag, comment); }

  /** Read particle positions and number of particles from XYZ file format.
   *  http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html
   *  [ # optional comment line ] comment line (can be blank)
      [ N                       ] # of atoms, required by this xyz reader
      [ molecule name           ] name of molecule (can be blank)
      atom1 x y z [optional data] atom name followed by xyz coords
      atom2 x y z [ ...         ] and and (optionally) other data. */
  void readXYZ(const char* fileName);

  // NOTE HWH: Depreciated for style. Use readXYZ.
  void readxyz(const char* fileName) { readXYZ(fileName); }

  /** Alternative readXYZ, this one adds and deletes molecules in addMolInit.
   *  See readxyz for description of xyz file format. */
  void readXYZ(std::ifstream& file);

  // NOTE HWH: Depreciated for style. Use readXYZ.
  void readxyz2(std::ifstream& file) { readXYZ(file); }

  /// Read particle possitions from XTC file format.
  #ifdef XDRFILE_H_
  int readXTC(const char* fileName, XDRFILE* trjFileXDR);

  /// Write particle possitions in XTC file format.
  void writeXTC(XDRFILE* trjFileXDR);
  #endif  // XDRFILE_H_

  /* Return the change in position according to periodic boundary conditions.
   * Assumes cubic periodic box centered about the origin.
   * NOTE HWH: Depreciated because it does not support triclinic box. */
  double pbc(const double x, const int dim);

  /** Return the change in position according to periodic boundary conditions.
   *  Assumes box centered about the origin, and that the particle only
   *  needs be wrapped once. */
  vector<double> pbc(const vector<double> x);

  /// Return random molecule as vector of particle numbers.
  vector<int> randMol();

  /** Return random molecule as vector of particle numbers.
   *  This is accomplished by picking a random atom, and finding all other
   *  atoms in the same molecule. In this case, the returned molecule must not
   *  be the same as the ones listed in jmpart. */
  vector<int> randMolDiff(const vector<int> jmpart);

  /// Single particle alternative to randMolDiff(vector).
  vector<int> randMolDiff(const int jMol);

  /** Return random molecule in jmpart as vector of particle numbers.
   *  This is accomplished by picking a random atom from jmpart, and finding
   *  all other atoms in the same molecule. */
  vector<int> randMolSubset(const vector<int> jmpart);

  /** Stores position of particles listed in mpart, and also orientations
   *  if not spherically symmetric. */
  void xStore(const vector<int> mpart);

  /// Store position and orientation of all particles.
  void xStoreAll();

  /** Store position of particles listed in mpart, and also orientations if
   *  not spherically symmetric. But for Multi implementation, store multiple
   *  instances of the coordinates before writing over them (e.g., for use
   *  with configurational bias). */
  void xStoreMulti(const vector<int> mpart, const int flag
   /** if flag == -1, clear all previous stores, then store particles in mpart
    * if flag == -2, store another mpart
    * if flag == positive integer, restore mpart particles from the 'index'th    *   store with negative flags. */
  );

  /** Restore list of particles, mpart, to the positions and orientations
   *  when the last xStore() was called. */
  void restore(const vector<int> mpart);

  /** Restore all particles to positions and orientations when last xStore()
   * was called. */
 void restoreAll();

  /// Set particle iPart to position "pos".
  void xset(double pos, int iPart, int dim) {x_[iPart*dimen_+dim] = pos; }

  /// Set particle iPart to position "pos".
  void xset(const int iPart, vector<double> pos)
    {for (int j = 0; j < static_cast<int>(pos.size()); ++j)
      x_[dimen_*iPart+j] = pos[j]; }

  /// Set length of domain boundary to "boxl" in given dimension.
  void lset(double boxl, int dimension) {l_[dimension] = boxl; }

  /// Set length of domain boundary to "boxl" in all dimensions.
  void lset(double boxl) {
    for (int dim = 0; dim < dimen_; ++dim) {l_[dim] = boxl; } }

  /** Randomly displace particles mpart by a random amount
   *  of maximum size maxDisp in each dimension. */
  void randDisp(const vector<int> mpart, const double maxDisp);

  // NOTE HWH: Depreciated function
  void randDisp(int part, const double maxDisp);

  // NOTE HWH: Depreciate function. Only used by TrialCluster, and no wrap?
  void randDispMulti(const vector<int> mpart, const double maxDisp);

  /** Random rotation of particles mpart by a random amount of maximum size
   *  maxRot.*/
  void randRotate(const vector<int> mpart, const double maxRot);

  /** Randomly rotate multiple molecules about the center of mass by a random
   *  amount given by maxDisp. Mass of all particles assumed equal
   *  to 1, or 0 if their "sig" size parameter is zero. */
  void randRotateMulti(const vector<int> mpart, const double maxDisp,
                       const vector<double> &sig);

  /// Alternative randRotateMulti where all particles have equal mass.
  void randRotateMulti(const vector<int> mpart, const double maxDisp);

  /// Translate molecule iMol by a displacement vector "r".
  void transMol(const int iMol, const vector<double> &r);

  /// Add particle with position v, type itype and molecule imol.
  void addPart(const vector<double> v, const int itype, const int imol);

  void addMolInit(const char* fileName);
  void addMolInit(const std::string fileName) { addMolInit(fileName.c_str()); }
  /// Add molecule of given type.
  void addMol(const char* type);

  /// Alternative addMol.
  void addMol(const std::string type) { addMol(type.c_str()); }

  /// Alternative addMol for index of order of addMolInits.
  void addMol(const int index = 0) { addMol(addMolListType_[index].c_str()); }

  /** Returns whether or not fast deletion method is applicable.
   *  Use a faster delete method if molecule of same type was the last
   *  molecule to be added.is at the by putting last molecule where mpart
   *  exists, and then delete the molecule at the end of the array.
   *  Also check particles to delete are part of entire molecule. */
  bool fastDelApplicable(const vector<int> mpart) const;

  /// Delete particle.
  void delPart(const int ipart);

  /// Delete list of particles.
  void delPart(const vector<int> mpart);

  /** Wrap molecule defined by list of particles, mpart, according to first
   *  particle in molecule. Assumes origin at the center of the boundary box.
   */
  void wrap(const vector<int> mpart);

  /** Wrap the position defined by rvec in simulation domain.
   *  Assumes origin at the center. */
  void rwrap(vector<double> *rvecPtr);

  /// Wrap all molecules.
  void wrapMol();

  /* Read bulk molecule xyz with natoms each, and assign type and mol.
   * NOTE HWH: Depreciated function. */
   void readXYZBulk(const int nMolAtoms, const char* type, const char* fileName);

  /// Returns particle IDs of the last molecule that was added.
  vector<int> lastMolIDVec();

  // Check the bond lengths of particles of type against reference particle.
  // NOTE HWH: Depreciated
  int checkBond(const char* type, const double tolerance);

  /** Check the bond lengths of particles of type against reference particle.
   *  \return 1 if bonds match. */
  int checkBond(const double tolerance);

  /// Generate atomic positions listed by molecule into variable xMol.
  void xMolGen();

  double minl() const;    //!< Returns minimum boundary length.

  /// Return whether atoms iAtom and jAtom are within the bond shell.
  int bonded(const int iAtom, const int jAtom, const double rAbove,
             const double rBelow);

  /** Alternative bonded where first particle in lists mpart and jmpart are
   *  applied to the bonded function described above. */
  int bonded(const vector<int> mpart, const vector<int> jmpart,
             const double rabove, const double rbelow);

  /** Moves molecule mpart to bonded/nonbonded region of jmpart
   *  and returns whether mpart was previously in a bonded configuration. */
  int avb(const vector<int> mpart, const vector<int> jmpart,
          const double rabove, const double rbelow,
          const char* type  //!< move to type=="(non)bonded" region
  );

  /// Move atom iAtom to in/out region of jmpart.
  // NOTE HWH : is this function still used?
  void avb(const int iAtom, const int jAtom, const double rAbove,
           const double rBelow, const char* region);

  /// Return squared distance between two points subject to PBCs.
  double rsq(const vector<double> xi, const vector<double> xj);

  /// Initialize cells in cell list
  void updateCells(const double dCellMin,  //!< minimm cell size
    const double rCut  //!< maximum cutoff radius for interactions
    );

  /// Simplified cell list initialization for dCellMin == rCut.
  void updateCells(const double dCellMin) { updateCells(dCellMin, dCellMin); }

  /** Further simplified cell list initialization for using previously stored
   *  value of dCellMin. */
  void updateCells() { updateCells(dCellMin_); }

  /// Assign particles or molecules to the cell list.
  void buildCellList();

  /// Return scalar cell index given particle number.
  int iatom2m(const double &ipart);

  /** Return scalar cell index given molecule number, which is based soley on
   *  the position of the first particle in the molecule. */
  int imol2m(const double &iMol) { return iatom2m(mol2part_[iMol]); }

  /// Generate neighbor list for iMol from cell list by molecule cutoff.
  void buildNeighListCell(const int iMol);

  /// Generate neighbor list for particle ipart from cell list by atom cuttoff.
  void buildNeighListCellAtomCut(const int ipart);

  void cellOff();                             //!< Turn off cell list.
  void updateCellofiMol(const int iMol);      //!< Updates cell for iMol.
  void updateCellofallMol();                  //!< Updates cell for all mols.

  /// Return 1 if no errors found in cell list.
  int checkCellList();

  /// Initialize cut-off method for cell list.
  // HWH NOTE: this is a bad name, because it often needs to be called even when
  // no cell list is involved
  void initCellAtomCut(const int flag);

  /* Position of origin to add next molecule called by addMol().
   * By default, xAdd is NULL which results in a random position. */
  vector<double> xAdd;

  /** Initialize quaternions of molecule iMol and stores current positions
   *  as reference positions. */
  void qMolInit(const int iMol);

  /** Initialize quaternions of all molecules and stores current positions
   *  as reference positions. */
  void qMolInit();

  /// Update current positions of molecule iMol using quaternions.
  void quat2pos(const int iMol);

  /// Update current positions of molecules imMol using quaternions.
  void quat2pos(const vector<int> imMol) {
    for (unsigned int i = 0; i < imMol.size(); ++i) quat2pos(imMol[i]);
  }

  /// Set particle type of iatom to itype.
  void settype(const int iatom, const int itype);

  /*  **  **  Do not include this in DOXYGEN  **  **  **
   * Initialize with LAMMPS data file. Read number of atoms, molecules,
   * atom types, masses. Resizes appropriate arrays. This function is
   * not be confused with addMolInit, which is preferred for users. */
  void initLMPData(const std::string fileName,
    const int nTypesExist = 0  // number of particle types that already exists
    );

  /* Initialize with data file, where file JSON file exensions are recognized
   * if JSON_ preprocessor macro is defined at comilation. */
  void initData(const std::string fileName, const int nTypesExist = 0);

  #ifdef JSON_
    // Initialize with JSON data file.
    void initJSONData(const std::string fileName, const int nTypesExist = 0);
  #endif  // JSON_

  /// Return 1 if size of protected arrays in this class pass all checks.
  int checkSizes();

  /// Add atom, iatom, as a tagged particle to track.
  void tagAtom(const int iatom) { tag_.push_back(iatom); }

  /// Remove the last atom from the tagged particles.
  void tagAtomPopBack() { tag_.pop_back(); tagStage_ = 0.; }

  /// Remove tag from all particles.
  void tagAtomClear() { tag_.clear(); tagStage_ = 0.; }

  /** Return list of all particles in molecule of the first tagged atom.
   *  Assumes only one tagged particle, which is first particle of mol. */
  vector<int> tag2mpart();

  /// Given addMolInit type, return space which contains only one molecule.
  shared_ptr<Space> findAddMolInList(const string typeStr);

  /// Given addMolInit type, return index of order in which addMolInit began.
  int findAddMolListIndex(const string typeStr);

  /// Change the bond lengths of molecule, iMol.
  void scaleMol(const int iMol, const vector<double> bondLengths);

  /// Change the bond lengths of molecule, iMol, while storing stage for tag.
  void scaleMol(const int iMol, const vector<double> bondLengths,
                const double stage)
                { tagStage_ = stage; scaleMol(iMol, bondLengths); }

  /// Update clusters of entire system.
  void updateClusters(const double rCut);

  /// Add particle type to consider in cluster analysis.
  void addTypeForCluster(const int type) { clusterType_.push_back(type); }

  /// Delete all particles of a given type.
  void delTypePart(const int type);

  /// Swap the particle coordinates of two objects with equal particle numbers.
  void swapPositions(Space *space);

  /** Swap positions of iMol and jMol. Currently only implemented for
   *  configurations with only monoatomic particles e.g., nMol == natom */
  void swapPositions(const int iMol, const int jMol);

  /// Return maximum distance between molecule center and atom in molecule.
  double maxMolDist();

  /// Return list of all particles in molecule iMol.
  vector<int> imol2mpart(const int iMol);

  /// Given list of particles, return inertia tensor (3D only).
  vector<vector<double> > inertialTensor(const vector<int> mpart);

  /// Return center of mass of list of particles. Assumes mass of 1.
  vector<double> rcom(const vector<int> mpart);

  /// Given list of particles, return list of molecules they comprise.
  vector<int> mpart2mmol(const vector<int> mpart);

  /// Print cluster statistics to file.
  void printClusterStat(const char* fileName);

  /** For each cluster, starting with first atom or molecule in cluster,
   *  find PBC of other cluster molecules and shift to generate xcluster(). */
  void xClusterGen();

  /// Compute shape metrics of clusters.
  void xClusterShape();

  /// Reset the stored cluster statistics.
  void clusterReset() { clusterSizeAccVec_.reset(); clusterNumAccVec_.reset();
                        clusterSizeDistribution_.reset(); }

  /// Update the stored cluster statistics.
  void updateClusterVars(const int nClusters);

  /// Use contact and contactpbc to update cluster variables.
  void contact2clusterAlt(
    /** For each molecule, list of molecules which are in contact. */
    vector<vector<int> > contact,
    /** For each molecule, PBC shift (for each dimension) of each contact. */
    vector<vector<vector<double> > > contactpbc
  );

  // Use contact and contactpbc to update cluster variables.
  // NOTE HWH: Deprciated. Use contact2clusterAlt.
  void contact2cluster(vector<vector<int> > contact,
                       vector<vector<vector<double> > > contactpbc);

  /// Place atom at the COM of all other atoms in list of particles, mpart.
  void setAtomAsCOM(const int atom, const vector<int> mpart);

  /// Place atom "iAtom" in sphere of radius "r" w.r.t. atom "jAtom"
  void setAtomInSphere(const int iAtom, const int jAtom, const double r);

  /** Place atom "iAtom" randomly in the circule of radius "r" w.r.t. atom
   *  "jAtom" and angle "theta" given by \f$ \angle ijk \f$. */
  void setAtomInCircle(const int iAtom, const int jAtom, const int kAtom,
                       const double r, const double theta);

  /** Place atom 3 in branch, given existing atoms a1, a2 and a4.
   *
   *               1
   *               .
   *       (t142)  .  (t143)
   *               4
   *             .   .
   *          .        .
   *        2    (t243) (3)
   */
  void setAtomInBranch(const int a1, const int a2, const int a3, const int a4,
    const double t143,  //!< angle between particles 1, 4 and 3
    const double t243,  //!< angle between particles 2, 4 and 3
    const double L      //!< bond length between particles 3 and 4
  );

  /** Modify bond angle of iAtom to theta, where bond angle is defined by
   *  angle \f$ \angle ijk \f$, preserving the plane that i,j,k reside. */
  void modBondAngle(const int iAtom, const int jAtom, const int kAtom,
                    const double theta);

  /// Modify all bond angles of type angleType in molType to angle theta.
  void modBondAngle(const int angleType, const double theta,
                    const char* molType);

  /** Return bond parameters (k,l0) for potential \f$ U=k(l-l0)^2 \f$ for
   *  atoms iAtom and jAtom. Returns 0 if non-existent. */
  vector<double> bondParams(const int iAtom, const int jAtom);

  /** Return angle parameters (k,t0) for potential \f$ U=k(t-t0)^2 \f$ for
   *  atoms \f$ \angle ijk \f$. Returns 0 if non-existent. */
  vector<double> angleParams(const int iAtom, const int jAtom, const int kAtom);

  /// Modify angle parameters.
  void modAngleParams(const int angleType, const int angleIndex,
    const double param) { angleParam_[angleType][angleIndex] = param; }

  /// Return list of bonds involving iAtom.
  vector<vector<int> > listBonds(const int iAtom);

  /// Return list of angles involving atoms iAtom and jAtom.
  vector<vector<int> > listAngles(const int iAtom, const int jAtom);

  // accumulate radial distance histogram
  // NOTE HWH: Depreciate in favor of Analysis class.
  void nRadialHist(Histogram *nhistPtr);

  // Write the radial distribution function.
  void printRadial(const Histogram &nhist, const char* fileName);

  /// Pivot iMol about the reflection point, r. \f$  rnew = 2*r - rAtom \f$.
  void pivotMol(const int iMol, const vector<double> r);

  /// Return a random position within the domain.
  vector<double> randPosition();

  /// Return a random position, centered about molecule iMol.
  vector<double> randPosition(const double iMol,
    const double maxDisp  //!< maximum distance away from iMol in each dimen
  );

  // compute the scattering intensity using full Debye equation
  // NOTE HWH: Depreciate in favor of Analysis.
  vector<double> scatterIntensity(const double qMin, const double qMax,
                                  const double dq);

  /** Scale the domain by a factor in dimension, dim. Scale positions of
   *  particles as well. If molecules are present, only scale the COM to
   *  maintain bonds lengths and angles. */
  void scaleDomain(const double factor, const int dim);

  /** Alternative scaleDomain which scales all dimensions by
   * \f$ factor^{1/dimen()} \f$, which is corresponds to scaling volume. */
  void scaleDomain(const double factor) { for (int dim = 0; dim < dimen_; ++dim)
    { scaleDomain(pow(factor, 1./dimen_), dim); } }

  /// Return the global, rotationally invariant q6 bond order parameter.
  double Q6(const double rCut  //!< distance cut-off to define neighbors
    );

  /// NOTE HWH: Define precisely XY, XZ and YZ tilt factors here.
  /// Set the xy tilt factor
  void setXYTilt(const double xyTilt);
  void setXZTilt(const double xzTilt);
  void setYZTilt(const double yzTilt);

  /// modify the xy tilt factor, and transform the particles
  // NOTE HWH CLEANUP: same for xy,yz,etc..
  void modXYTilt(const double deltXYTilt);
  void modXZTilt(const double deltXZTilt);
  void modYZTilt(const double deltYZTilt);

  /// Return minimum bond length in all molecules present or in addMolInit.
  double minBondLength();

  /// Initialize euler angle representation for orientation.
  void initEuler(const int flag) { eulerFlag_ = flag; }

  /// update the euler angles of iMol according to position relative to ref
//  void pos2euler(const int iMol);

  /// Return the index of a randomly selected moleule of type iMolType.
  int randMolofType(const int iMolType);

  /// Set maximum box length.
  void setMaxBoxLength() { maxlFlag_ = 1; maxl_ = l_; }

  /// Print vmd script for visualization of trajectories.
  void printxyzvmd(const char* fileName,
    /** If initFlag == 1, format for XYZ files.
     *  If initFlag == 2, format for XTC files. */
    const int initFlag);

  /// Initialize equimolar constraint.
  void equiMolar(
    /** If flag == 1, double branch//add or delete any on even nMol.
     *  If flag == 2, single branch//start with iMolType 0.
     *  If flag == 3, single branch//start with iMolType 1. */
    const int flag) { equiMolar_ = flag; }

  /** Given particle ipart, return euler angles. This routine assumes that
   *  the anisotropic particle is a solid of revolution, and the axis of
   *  symmetry points along the ipart->ipart+1 unit vector. */
  vector<double> ipart2euler(const int ipart);

  /// replicate the system via peroidic boundary conditions
  void replicate(
    const int nx = 1,  //!< number of times to replicate in x dimension
    const int ny = 1,  //!< number of times to replicate in y dimension
    const int nz = 1   //!< number of times to replicate in z dimension
    );

  /// Initialize intramolecular interactions via contact map (1 if ixn).
  void initIntra(const vector<vector<int> >& map);

  // functions for read-only access of private data-members
  /// full access to private data-members
  vector<vector<vector<int> > > intraMap() { return intraMap_; }
  int dimen() const { return dimen_; }
  int qdim() const { return qdim_; }
  int id() const { return id_; }
  int natom() const { return static_cast<int>(x_.size())/dimen_; }
  vector<double> x() const { return x_; }
  vector<double> xcluster() const { return xcluster_; }
  vector<vector<vector<double> > > xMol() const { return xMol_; }
  vector<vector<double> > xold() const { return xold_; }
  vector<vector<vector<double> > > xOldMulti() const { return xOldMulti_; }
  double x(int ipart, int dim) const { return x_[dimen_*ipart+dim]; }
  vector<int> mol2part() const { return mol2part_; }
  vector<int> tag() const { return tag_; }
  double tagStage() const { return tagStage_; }
  vector<double> l() const { return l_; }
  double l(const int i) const { return l_[i]; }
  double vol() const { return product(l_); }  //!< simulation domain volume
  double type(const int i) const { return type_[i]; }
  vector<int> type() const { return type_; }
  vector<int> mol() const { return mol_; }
  vector<int> nType() const { return nType_; }
  vector<int> nMolType() const { return nMolType_; }
  vector<int> listAtoms() const { return listAtoms_; }
  vector<int> listMols() const { return listMols_; }
  int nParticleTypes() const { return static_cast<int>(nType_.size()); }
  int nMolTypes() const { return static_cast<int>(nMolType_.size()); }
  vector<int> nCellVec() const { return nCellVec_; }
  int nCell() const { return nCell_; }
  vector<string> moltype() const { return moltype_; }
  vector<int> molid() const { return molid_; }
  vector<vector<vector<double> > > xMolRef() const { return xMolRef_; }
  vector<double> qMol() const { return qMol_; }
  double qMol(const int iMol, const int dim) const
    { return qMol_[qdim_*iMol+dim]; }
  vector<double> qMol(const int iMol) const;
  void qMolAlt(const int iMol, const int dim, const double q)
    { qMol_[qdim_*iMol+dim] = q; }
  bool fastDel() const { return fastDel_; }
  int fastDelMol() const { return fastDelMol_; }
  int cellType() const { return cellType_; }
  vector<vector<int> > neighCell() const { return neighCell_; }
  vector<int> neighListCell() const { return neighListCell_; }
  vector<int> neighListChosen() const { return *neighListChosen_; }
  vector<vector<int> > cellList() const { return cellList_; }
  vector<int> atom2cell() const { return atom2cell_; }
  double dCellMin() const { return dCellMin_; }
  int nMol() const { return static_cast<int>(moltype_.size()); }
  vector<vector<int> > cMaskPnt() const { return cMaskPnt_; }
  vector<shared_ptr<Space> > addMolList() const { return addMolList_; }
  vector<string> addMolListType() const { return addMolListType_; }
  shared_ptr<std::ifstream> xyzFile() const { return xyzFile_; }
  bool xyzFileEOF() const
    { if (xyzFile_ == NULL) { return false; } else { return xyzFile_->eof();} }
  vector<int> cluster() const { return cluster_; }
  vector<int> clusterMol() const { return clusterMol_; }
  vector<int> clusterSizes() { return clusterSizes_; }
  int nClusters() const { return static_cast<int>(clusterSizes_.size()); }
  double clusterAvSize() const {
    if (static_cast<int>(clusterSizes_.size()) == 0) { return 0;
    } else { return vecAverage(clusterSizes_); } }
  vector<int> clusterType() const { return clusterType_; }
  vector<vector<int> > clusterList() const { return clusterList_; }
  AccumulatorVec clusterSizeAccVec() const { return clusterSizeAccVec_;}
  AccumulatorVec clusterNumAccVec() const { return clusterNumAccVec_;}
  AccumulatorVec clusterSizeDistribution() const
    { return clusterSizeDistribution_;}
  Accumulator freeMon() const { return freeMon_; }
  double clusterAsphericityAv() const {
    if (static_cast<int>(clusterAsphericity_.size()) == 0) { return 0;
    } else { return vecAverage(clusterAsphericity_); } }
  double clusterAcylindricityAv() const {
    if (static_cast<int>(clusterAcylindricity_.size()) == 0) { return 0;
    } else { return vecAverage(clusterAcylindricity_); } }
  double clusterRelShapeAnisoAv() const {
    if (static_cast<int>(clusterRelShapeAniso_.size()) == 0) { return 0;
    } else { return vecAverage(clusterRelShapeAniso_); } }
  double clusterRgAv() const {
    if (static_cast<int>(clusterRg_.size()) == 0) { return 0;
    } else { return vecAverage(clusterRg_); } }
  void preMicellarAgg(const int size) { preMicellarAgg_ = size; }
  vector<vector<int> > bondList() const { return bondList_; }
  vector<vector<double> > bondParam() const { return bondParam_; }
  vector<vector<int> > angleList() const { return angleList_; }
  vector<vector<double> > angleParam() const { return angleParam_; }
  bool sphereSymMol() const { return sphereSymMol_; }
  double xyTilt() const { return xyTilt_; }
  double xzTilt() const { return xzTilt_; }
  double yzTilt() const { return yzTilt_; }
  int floppyBox() const { return floppyBox_; }
  int equiMolar() const { return equiMolar_; }
  int percolation() const { return percolation_; }
  int eulerFlag() const { return eulerFlag_; }

 private:
  int dimen_;     //!< dimesion of real space
  int id_;   //!< ID of space class
  int qdim_;     //!< dimesion of quaternion space (dimen_ + 1)
  vector<double> x_;        //!< atomic positions

  /// Set the default values during construction.
  void defaultConstruction_();

  /// Reconstruct pointers pointers upon cloning.
  void reconstruct_();

  /// returns first atom for each molecule. sorted same as x_. number of
  //  elements is nMol+1 where last element is natom()
  vector<int> mol2part_;

  /// atomic positions by molecule[mol][partInMol][dim]
  vector<vector<vector<double> > > xMol_;
  vector<int> type_;        //!< atomic type

  /// for a given atomic type, return total number of particles with that type
  vector<int> nType_;

  /// for a given molecule type, return total number of molecules with that type
  vector<int> nMolType_;
  vector<int> mol_;     //!< atomic molecular id
  vector<string> moltype_;  //!< moltype_[mol] = string type of molecule
  vector<int> molid_;   //!< molid_[mol] = integer type of molecule
  vector<int> tag_;     //!< tag atom index, update with insertions or deletions
  double tagStage_;     //!< stage of tagged atom
  vector<vector<double> > xold_;  //!< old atomic positions before randDisp
  vector<double> xOldAll_;        //!< old atomic positions

  /// stores series of old atomic positions
  vector<vector<vector<double> > > xOldMulti_;
  vector<int> listAtoms_;  //!< list of consequetive integers (0 to nAtom() - 1)
  vector<int> listMols_;   //!< list of consequetive integers (0 to nMol() - 1)
  vector<int> cluster_;    //!< cluster id of particles (not auto updated)
  vector<int> clusterMol_;    //!< cluster id of molecules (not auto updated)
  vector<int> clusterSizes_;  //!< number of particles in each cluster id
  vector<double> xcluster_;   //!< position of clusters not wrapped in pbc

  /// list of atom types to consider in clustering algorithm
  vector<int> clusterType_;

  /// for each cluster type, list of particles
  vector<vector<int> > clusterList_;
  AccumulatorVec clusterSizeAccVec_;    //!< accumulator for cluster sizes
  AccumulatorVec clusterNumAccVec_;     //!< accumulator for number of clusters

  /// accumulator for cluster size dist
  AccumulatorVec clusterSizeDistribution_;
  Accumulator freeMon_;            //!< pre-micellar aggregate concentration
  vector<double> clusterAsphericity_;    //!< asphericity of each cluster
  vector<double> clusterAcylindricity_;  //!< acylindricity of each cluster
  vector<double> clusterRelShapeAniso_;  //!< relative shape anisotropy
  vector<double> clusterRg_;             //!< radius of gyration of each cluster
  int preMicellarAgg_;   //!< cluster size as cut-off for premicellar aggregates
  int percolation_;      //!< flag if percolation was detected

  /// flood fill algorithm to identify clusters based on atomic distance cutoff
  void floodFill3d_(const int clusterNode, const int clusterID,
                   const double rCut);
  void floodFillCell3d_(const int clusterNode, const int clusterID,
                       const double rCut);
  void floodFill2d_(const int clusterNode, const int clusterID,
                   const double rCut);
  void floodFillContact_(const int clusterNode, const int clusterID,
                        vector<vector<int> > *contactPtr,
                        vector<vector<vector<double> > > *contactpbcPtr);
  void floodFillContactAlt_(const int clusterNode, const int clusterID,
                           vector<vector<int> > *contactPtr,
                           vector<vector<vector<double> > > *contactpbcPtr,
                           vector<vector<int> > *image);

  /// prefile cluster vars
  void prefilClusterVars_();

  bool fastDel_;         //!< use fast method of deleting particles
  int fastDelMol_;       //!< molecule last deleted by fast method

  /// use cell list of type cellType,
  //   none = 0, AllenTildesley = 1, Hunenberger = 2
  int cellType_;
  bool cellAtomCut_;     //!< use atom cut or cell cut
  double dCellMin_;     //!< minimum cell size
  vector<int> nCellVec_;  //!< number of cells in each dimension
  int nCell_;               //!< number of cells
  vector<double> dCell_;  //!< width of cells in each dimension
  /// given cell id, obtain beginning and end cells for each stripe
  vector<vector<int> > cMaskPnt_;
  vector<vector<int> > cellList_;   //!< for cell, lists molecules in cell
  vector<int> mol2cell_;   //!< for given molecule, list cell
  vector<int> atom2cell_;   //!< for given atom, list cell
  /// for a given icell, neighbors of icell = neighCell_[icell]
  vector<vector<int> > neighCell_;
  vector<int> neighListCell_;       //!< generated neighbor list from cell list
  /// choosen neighlist (from neighListCell or all (listMols)
  vector<int> *neighListChosen_;

  /// Return the corners of cell "m".
  vector<vector<double> > cellCorners_(const int m);

  /// Return grid position of cell given cell index.
  vector<int> m2vec_(const int cellIndex);

  /// Return scalar index of cell given spatial coordinates.
  int rvec2m_(const vector<double> &r);

  /// Return scalar cell index given 3D grid coordinates (i, j, k).
  int mvec2m3d_(const int &i, const int &j, const int &k) const;

  /// Return scalar cell index given 2D grid coordinates (i, j).
  int mvec2m2d_(const int &i, const int &j) const;

  void eraseMolFromCell_(const int iMol);     //!< removes molecule from cell
  void addMoltoCell_(const int iMol);         //!< adds molecule to cellList_
  void eraseAtomFromCell_(const int ipart);   //!< removes atom from cellList_
  void addAtomtoCell_(const int ipart);       //!< adds atom to cellList_


  /// are molecules spherically symmetric, no rotations or quaternions necessary
  bool sphereSymMol_;
  vector<double> qMol_;    //!< orientation of molecules via quaternions
  /// old orientation of molecules via quaternions
  vector<vector<double> > qMolOld_;
  vector<double> qMolOldAll_;  //!< old orientation of molecules via quaternions
  /// multiple old orientation of molecules via quaternions
  vector<vector<vector<double> > > qMolOldMulti_;
  /// reference vector of molecules, xMolRef[mol][atom][dim]
  vector<vector<vector<double> > > xMolRef_;
  /// old reference vector of molecules
  vector<vector<vector<double> > > xMolRefOld_;
  int eulerFlag_;          //!< flag to use euler angles instead of quaternions
  /// list of molecules that may be added to simulation
  vector<shared_ptr<Space> > addMolList_;
  /// type of molecule that is listed in addMolList
  vector<string> addMolListType_;

  // i-o
  /// pointer xyz file to keep open while reading
  shared_ptr<std::ifstream> xyzFile_;

  // intramolecular interactions
  vector<vector<double> > bondParam_;  //!< bond parameters for each bond type
  vector<vector<double> > angleParam_;  //!< angle parameters for each angletype
  /// list of bonds[bond#] = [type, atom1, atom2]
  vector<vector<int> > bondList_;
  /// list of angles[angle#] = [type, a1, a2, a3]
  vector<vector<int> > angleList_;

  // variables which describe the domain
  /// simulation domain length for real space periodic boundaries
  vector<double> l_;
  /// flag is 1 if domain is constant throughout simulation, 0 otherwise
  double xyTilt_;                  //!< xy tilt factor
  double xzTilt_;                  //!< xz tilt factor
  double yzTilt_;                  //!< yz tilt factor
  /// flag is 1 if the box tilt factors are changing during the simulation
  int floppyBox_;
  vector<double> maxl_;     //!< max box length
  int maxlFlag_;            //!< flag if maxl set

  int equiMolar_;    //!< flag to enforce equimolar

  /// intraMap is a matrix of atom numbers that are 1 if they interact
  //  for multiple molecule types, addmol list in space can be used
  //  to call the map for individual molecule types
  vector<vector<vector<int> > > intraMap_;

  /// solve equations for branch, returning particle 3 given 1 and 2
  void solveBranch_(const double x1, const double y1, const double z1,
                   const double x2, const double y2, const double z2,
                   double *x3, double *y3, double *z3, const double c143,
                   const double c243);

};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // SRC_SPACE_H_
