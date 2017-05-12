/**
 * \brief coordinates
 *
 * The space class owns variables and functions associated with the real-space
 * position of atoms and the domain.
 * The domain is a cubic box centered about the origin, subject to periodic
 * boundary conditions
 */

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

namespace feasst {

class Space : public BaseAll {
 public:
  Space(int dimen = 3, int id = 0);
  explicit Space(const char* fileName);
  ~Space();
  Space* clone() const;
  shared_ptr<Space> cloneShrPtr() const;

  // defaults in constructor
  void defaultConstruction();

  // reconstruct pointers pointers upon cloning
  void reconstruct();

  /// write restart file
  void writeRestart(const char* fileName);

  //// initialize configuration with natom atoms
  int init_config(const int natom);

  /// write configuration
  int printxyz(const char* fileName, const int initFlag,
    const std::string comment="");

  /// read particle positions and number of particles from XYZ file format
  void readxyz(const char* fileName);
  void readxyz2(std::ifstream& file);  // alternative method used addMolInit

  /// read particle possitions from XTC file format
  #ifdef XDRFILE_H_
  int readXTC(const char* fileName, XDRFILE* trjFileXDR);
  void writeXTC(XDRFILE* trjFileXDR);
  // int readXTC(const char* fileName);
  #endif  // XDRFILE_H_

  /// periodic boundary conditions. For orthologal box vectors only
  double pbc(const double x, const int i);

  /// periodic boundary conditions
  vector<double> pbc(const vector<double> x);

  /// return random molecule as vector of particle numbers, mpart
  vector<int> randMol();

  //// randomly select molecule not including particles jmpart
  vector<int> randMolDiff(const vector<int> jmpart);
  vector<int> randMolDiff(const int jMol);  //!< not including molecule jMol

  /// randomly select only from particles in jmpart
  vector<int> randMolSubset(const vector<int> jmpart);

  /// restore configuration to original positions, stored in xold_ for mpart
  void restore(const vector<int> mpart);
  void restoreAll();

  /// set particle positions
  void xset(double c, int iAtom, int dim) {x_[iAtom*dimen_+dim] = c; }
  void xset(const int ipart, vector<double> v)
    {for (int j = 0; j < static_cast<int>(v.size()); ++j)
      x_[dimen_*ipart+j] = v[j]; }

  /// set domain boundaries
  void lset(double c, int i) {l_[i] = c; }
  void lset(double c) { for (int dim = 0; dim < dimen_; ++dim) {l_[dim] = c; } }

  /// random displacement of particle(s)
  void randDisp(const vector<int> mpart, const double maxDisp);
  void randDispMulti(const vector<int> mpart, const double maxDisp);
  void randRotate(const vector<int> mpart, const double maxDisp);
  void randRotateMulti(const vector<int> mpart, const double maxDisp);
  void randRotateMulti(const vector<int> mpart, const double maxDisp,
                       const vector<double> &sig);
  void randDisp(int part, const double maxDisp);

  /// translate molecule
  void transMol(const int iMol, const vector<double> &r);

  /// delete particle
  void delPart(const int ipart);
  void delPart(const vector<int> mpart);

  /// add one particle
  void addPart(const vector<double> v, const int itype, const int imol);

  /// add molecule of predefined type
  void addMol(const char* type);
  void addMol(const std::string type) { addMol(type.c_str()); }
  void addMolInit(const char* fileName);
  void addMolInit(const std::string fileName) { addMolInit(fileName.c_str()); }
  void addMol(const int index) { addMol(addMolListType_[index].c_str()); }
  void addMol() { addMol(0); }

  // ensure particles are within domain boundaries
  /// wrap molecule defined by mpart according to first particle in molecule
  void wrap(const vector<int> mpart);

  /// wrap position defined by rvec in simulation domain
  void rwrap(vector<double> *rvecPtr);
  void wrapMol();     //!< wrap all molecules

  /// read bulk molecule xyz with natoms each, and assign type and mol
  void readXYZBulk(const int nMolAtoms, const char* type, const char* fileName);

  /// returns vector of particle IDs of last molecule
  vector<int> lastMolIDVec();

  /// checks bond lengths of particles against reference particle xRef
  int checkBond(const char* type, const double tol);
  int checkBond(const double tol);

  /// generate atomic positions listed by molecule into variable xMol_
  void xMolGen();

  double minl() const;    //!< returns minium boundary distance

  /// returns whether molecules mpart and jmpart are bonded
  int bonded(const vector<int> mpart, const vector<int> jmpart,
             const double rabove, const double rbelow);
  int bonded(const int iAtom, const int jAtom, const double rAbove,
             const double rBelow);

  /// moves molecule mpart to bonded/nonbonded region of jmpart,
  //  and returns whether mpart was previously in a bonded configuration
  int avb(const vector<int> mpart, const vector<int> jmpart,
          const double rabove, const double rbelow, const char* type);

  /// moves atom iAtom to in/out region of jAtom
  void avb(const int iAtom, const int jAtom, const double rAbove,
           const double rBelow, const char* region);

  /// stores position of a vector of particles in xold_
  void xStore(const vector<int> mpart);
  void xStoreAll();

  /// stores position of a vector of particles in xOldMulti_
  void xStoreMulti(const vector<int> mpart, const int flag);

  /// return squared distance between two points subject to pbc
  double rsq(const vector<double> xi, const vector<double> xj);

  /// returns whether or not fast deletion method is applicable
  bool fastDelApplicable(const vector<int> mpart) const;

  /// cell list
  void updateCells(const double dCellMin, const double rCut);
  void updateCells(const double dCellMin) { updateCells(dCellMin, dCellMin); }
  void updateCells() { updateCells(dCellMin_); }
  void buildCellList();
  vector<vector<double> > cellCorners(const int m);  //!< return corners of cell
  vector<int> m2vec(const int m);  //!< return vector position of cell

  /// return scalar cell index given position
  int rvec2m(const vector<double> &r);

  /// return scalar cell index given vector cell index
  int mvec2m3d(const int &i, const int &j, const int &k) const;

  /// return scalar cell index given vector cell index
  int mvec2m2d(const int &i, const int &j) const;

  /// return scalar cell index given molecule number
  int imol2m(const double &iMol);

  /// return scalar cell index given atom number
  int iatom2m(const double &ipart);

  /// generate neighbor list for iMol from cell list
  void buildNeighListCell(const int iMol);

  /// generate neighbor list for ipart from cell list
  void buildNeighListCellAtomCut(const int ipart);
  void cellOff();                             //!< turn off cell list
  void eraseMolFromCell(const int iMol);      //!< removes molecule from cell
  void addMoltoCell(const int iMol);          //!< adds molecule to cellList_
  void eraseAtomFromCell(const int ipart);    //!< removes atom from cellList_
  void addAtomtoCell(const int ipart);        //!< adds atom to cellList_
  void updateCellofiMol(const int iMol);      //!< updates cell for iMol
  void updateCellofallMol();                  //!< updates cell for all mols

  /// returns 1 if no errors found in cell list. stores current cell list,
  //  rebuilds, and compares
  int checkCellList();

  /// initialize cut-off method for cell list
  // HWH NOTE: this is a bad name, because it often needs to be called even when no
  // cell list is involved
  void initCellAtomCut(const int flag);

  /// if !null, position to add center of molecule in addMol
  vector<double> xAdd;

  /// initialize quaternions of all molecules
  //  and stores current positions as reference positions
  void qMolInit();
  void qMolInit(const int iMol);

  /// update current positions using quaternions, given molecule numbers
  void quat2pos(const vector<int> imMol);
  void quat2pos(const int iMol);

  /// set particle type
  void settype(const int iatom, const int itype);

  /// initialize with data file
  void initData(const std::string fileName, const int nTypesExist = 0);

  /// initialize with LAMMPS data file
  void initLMPData(const std::string fileName, const int nTypesExist = 0);

  /// initialize with JSON data file
  #ifdef JSON_
    void initJSONData(const std::string fileName, const int nTypesExist = 0);
  #endif  // JSON_

  /// checks array sizes
  int checkSizes();

  /// tag atom, and update its index with insertions, deletions, or sorting
  void tagAtom(const int iatom);
  void tagAtomPopBack();
  void tagAtomClear() { tag_.clear(); tagStage_ = 0.; }
  vector<int> tag2mpart();

  /// find pointer to space of addMol in addMolList
  shared_ptr<Space> findAddMolInList(const string typeStr);
  int findAddMolListIndex(const string typeStr);

  /// scale molecule
  void scaleMol(const int iMol, const vector<double> bondLengths);
  void scaleMol(const int iMol, const vector<double> bondLengths,
                const double stage)
                { tagStage_ = stage; scaleMol(iMol, bondLengths); }

  /// flood fill algorithm to identify clusters based on atomic distance cutoff
  void floodFill3d(const int clusterNode, const int clusterID,
                   const double rCut);
  void floodFillCell3d(const int clusterNode, const int clusterID,
                       const double rCut);
  void floodFill2d(const int clusterNode, const int clusterID,
                   const double rCut);
  void floodFillContact(const int clusterNode, const int clusterID,
                        vector<vector<int> > *contactPtr,
                        vector<vector<vector<double> > > *contactpbcPtr);
  void floodFillContactAlt(const int clusterNode, const int clusterID,
                           vector<vector<int> > *contactPtr,
                           vector<vector<vector<double> > > *contactpbcPtr,
                           vector<vector<int> > *image);

  /// update clusters of entire system
  void updateClusters(const double rCut);

  /// add type for cluster analysis
  void addTypeForCluster(const int type) { clusterType_.push_back(type); }

  /// delete all particles of a given type
  void delTypePart(const int type);

  /// swap positions
  void swapPositions(Space *space);
  void swapPositions(const int iMol, const int jMol);

  /// maximum distance between molecule center and atom in molecule
  double maxMolDist();

  /// given molecule number, return vector of particles in molecule
  vector<int> imol2mpart(const int iMol);

  /// given list of particles, compute inertia tensor
  vector<vector<double> > inertialTensor(const vector<int> mpart);

  /// given list of particles, compute center of mass
  vector<double> rcom(const vector<int> mpart);

  /// given list of particles, obtain list of molecules
  vector<int> mpart2mmol(const vector<int> mpart);

  /// print cluster statistics
  void printClusterStat(const char* fileName);

  /// generate xcluster
  void xClusterGen();

  /// generate shape metrics of clusters
  void xClusterShape();

  /// reset cluster statistics
  void clusterReset() { clusterSizeAccVec_.reset(); clusterNumAccVec_.reset();
                        clusterSizeDistribution_.reset(); }

  /// update Cluster Vars
  void updateClusterVars(const int nClusters);

  /// prefile cluster vars
  void prefilClusterVars();

  /// use contact and contactpbc to update cluster variables
  void contact2cluster(vector<vector<int> > contact,
                       vector<vector<vector<double> > > contactpbc);
  void contact2clusterAlt(vector<vector<int> > contact,
                          vector<vector<vector<double> > > contactpbc);

  /// place atom at the COM of all other atoms in mpart
  void setAtomAsCOM(const int atom, const vector<int> mpart);

  /// place atom i in sphere of radius r w.r.t. atom j
  void setAtomInSphere(const int iAtom, const int jAtom, const double r);

  /// place atom i in circule of radius r w.r.t. atom j and angle theta
  //  w.r.t atom j and k
  void setAtomInCircle(const int iAtom, const int jAtom, const int kAtom,
                       const double r, const double theta);

  /// place atom 3 in branch, given theta143, theta243, and bond length l
  void setAtomInBranch(const int a1, const int a2, const int a3, const int a4,
                       const double t143, const double t243, const double L);

  /// modify bond angle of iAtom to theta, where bond angle is defied by
  //  angle <ijk, preserving the plane that i,j,k reside
  void modBondAngle(const int iAtom, const int jAtom, const int kAtom,
                    const double theta);
  void modBondAngle(const int angleType, const double theta,
                    const char* molType);

  /// solve equations for branch, returning particle 3 given 1 and 2
  void solveBranch(const double x1, const double y1, const double z1,
                   const double x2, const double y2, const double z2,
                   double *x3, double *y3, double *z3, const double c143,
                   const double c243);

  /// return bond parameters (U=k(l-l0)^2) for atoms i and j (0 if non-existant)
  vector<double> bondParams(const int iAtom, const int jAtom);

  /// return angle parameters (U=k(l-t0)^2) for atoms i,j,k (0 if non-existant)
  vector<double> angleParams(const int iAtom, const int jAtom, const int kAtom);

  /// modify angle parameters
  void modAngleParams(const int angleType, const int angleIndex,
    const double param) { angleParam_[angleType][angleIndex] = param; }

  /// return list of bonds involving iAtom
  vector<vector<int> > listBonds(const int iAtom);

  /// return list of angles involving atoms i and j
  vector<vector<int> > listAngles(const int iAtom, const int jAtom);

  /// accumulate radial distance histogram
  void nRadialHist(Histogram *nhistPtr);

  /// write the radial distribution function
  void printRadial(const Histogram &nhist, const char* fileName);

  /// pivot iMol about the reflection point, r
  void pivotMol(const int iMol, const vector<double> r);

  /// return a random position within the domain
  vector<double> randPosition();
  vector<double> randPosition(const double iMol, const double maxDisp);

  /// compute the scattering intensity using full Debye equation
  vector<double> scatterIntensity(const double qMin, const double qMax,
                                  const double dq);

  /// scale the domain by a factor
  void scaleDomain(const double factor) { for (int dim = 0; dim < dimen_; ++dim)
    { scaleDomain(pow(factor, 1./dimen_), dim); } }
  void scaleDomain(const double factor, const int dim);

  /// compute the global, rotationally invariant q6 bond order parameter
  double Q6(const double rCut);

  /// set the xy tilt factor
  void setXYTilt(const double xyTilt);
  void setXZTilt(const double xzTilt);
  void setYZTilt(const double yzTilt);

  /// modify the xy tilt factor, and transform the particles
  void modXYTilt(const double deltXYTilt);
  void modXZTilt(const double deltXZTilt);
  void modYZTilt(const double deltYZTilt);

  /// return minimum bond length in molecule
  double minBondLength();

  /// initialize euler angle representation for orientation
  void initEuler(const int flag) { eulerFlag_ = flag; }
  int eulerFlag() const { return eulerFlag_; }

  /// update the euler angles of iMol according to position relative to ref
//  void pos2euler(const int iMol);

  /// randomly select a molecule of a given type
  int randMolofType(const int iType);

  /// set maximum box length
  void setMaxBoxLength() { maxlFlag_ = 1; maxl_ = l_; }

  /// print vmd script for xyz files
  void printxyzvmd(const char* fileName, const int initFlag);

  /// initialize equimolar constraint
  // if flag == 1, double branch//add or delete any on even nMol
  // if flag == 2, single branch//start with iMolType 0
  // if flag == 3, single branch//start with iMolType 1
  void equiMolar(const int flag) { equiMolar_ = flag; }

  /// given ipart, compute euler angle
  vector<double> ipart2euler(const int ipart);

  /// replicate the system via peroidic boundary conditions
  void replicate(
    const int nx = 1,  //!< number of times to replicate in x dimension
    const int ny = 1,  //!< number of times to replicate in y dimension
    const int nz = 1   //!< number of times to replicate in z dimension
    );
  
  /// full access to private data-members
  vector<vector<vector<int> > > intraMap() { return intraMap_; }
  void initIntra(const vector<vector<int> >& map);

  // functions for read-only access of private data-members
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
  double vol() const { return feasst::product(l_); }  //!< simulation domain volume
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
  int nMol() const { return static_cast<int>(mol2part_.size() - 1); }
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
  
 private:
  int dimen_;     //!< dimesion of real space
  int id_;   //!< ID of space class
  int qdim_;     //!< dimesion of quaternion space (dimen_ + 1)
  vector<double> x_;        //!< atomic positions

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
};

}  // namespace feasst

#endif  // SRC_SPACE_H_
