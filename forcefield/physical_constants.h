#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/// constants from CODATA
/// http://dx.doi.org/10.1103/RevModPhys.84.1527
/// http://dx.doi.org/10.1063/1.4724320
const double boltzmannConstant = 1.3806488E-23; //J/K
const double avogadroConstant = 6.02214129E+23; //1/mol
const double elementaryCharge = 1.602176565E-19;  //C
const double permitivityVacuum = 8.854187817E-12;  //C^2/J/m
const double idealGasConstant = boltzmannConstant*avogadroConstant;
const double joulesPercal = 4.184;

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
