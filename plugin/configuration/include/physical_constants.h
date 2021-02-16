
#ifndef FEASST_CONFIGURATION_PHYSICAL_CONSTANTS_H_
#define FEASST_CONFIGURATION_PHYSICAL_CONSTANTS_H_

#include <sstream>
#include <memory>
#include <map>
#include <string>
#include "utils/include/arguments.h"

namespace feasst {

class PhysicalConstants {
 public:
  PhysicalConstants() {}

  //@{
  /** @name Base constants
   */

  /// Boltzman constant in units of Joules per Kelvin
  virtual const double boltzmann_constant() const = 0;

  /// Avogadro's constant units of number of particles per mol
  virtual const double avogadro_constant() const = 0;

  /// Permitivity of vacuum in units of C^2/J/m
  virtual const double permitivity_vacuum() const = 0;

  /// Elementary charge in units of C
  virtual const double elementary_charge() const = 0;

  //@}
  /** @name Derivatives
    Constants derived from the base constants.
   */
  //@{

  /// Ideal gas constant in units of Joules per Kelvin per mol
  const double ideal_gas_constant() const { return ideal_gas_constant_; }

  /// Convert e^2/Angstrom to kJ/mol by a factor of units (kJ*A/e^2/mol)
  const double charge_conversion() const { return charge_conversion_; }

  //@}
  /** @name Conversions
    Utility functions for unit conversion.
   */
  //@{

  /// Convert the number density (1/A^3) to g/cm^3
  inline double number_density_to_grams_per_cm3(const double density, /// 1/A^3
      const double molecular_weight) {
    return density*molecular_weight/avogadro_constant()*1e24; }

  /// Convert the density in g/cm^3 to number density (1/A^3)
  inline double grams_per_cm3_to_number_density(const double density, /// g/cm^3
      const double molecular_weight) {
    return density/number_density_to_grams_per_cm3(1., molecular_weight); }

  //@}

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const = 0;
  virtual std::shared_ptr<PhysicalConstants> create(
    std::istream& istr) const = 0;
  std::map<std::string, std::shared_ptr<PhysicalConstants> >& deserialize_map();
  std::shared_ptr<PhysicalConstants> deserialize(std::istream& istr);
  explicit PhysicalConstants(std::istream& istr) {
    istr >> class_name_; }
  virtual ~PhysicalConstants() {}

 protected:
  std::string class_name_ = "PhysicalConstants";
  void compute_derived_();

 private:
  double ideal_gas_constant_;
  double charge_conversion_;
};

/**
 * Physical constants from CODATA 2018
 */
class CODATA2018 : public PhysicalConstants {
 public:
  CODATA2018() : PhysicalConstants() {
    class_name_ = "CODATA2018";
    compute_derived_(); }

  const double boltzmann_constant() const override { return 1.380649E-23; }
  const double avogadro_constant() const override { return 6.02214076E+23; }
  const double permitivity_vacuum() const override { return 8.8541878128E-12; }
  const double elementary_charge() const override { return 1.602176634E-19; }

  std::shared_ptr<PhysicalConstants> create(std::istream& istr) const override {
    return std::make_shared<CODATA2018>(istr); }
  void serialize(std::ostream& ostr) const override {
    ostr << class_name() << " "; }
  explicit CODATA2018(std::istream& istr);
  virtual ~CODATA2018() {}
};

inline std::shared_ptr<CODATA2018> MakeCODATA2018() {
  return std::make_shared<CODATA2018>();
}

/**
 * Physical constants from CODATA 2014
 * http://dx.doi.org/10.1103/RevModPhys.88.035009
 * http://dx.doi.org/10.1063/1.4954402
 */
class CODATA2014 : public PhysicalConstants {
 public:
  CODATA2014() : PhysicalConstants() {
    class_name_ = "CODATA2014";
    compute_derived_(); }

  const double boltzmann_constant() const override { return 1.38064852E-23; }
  const double avogadro_constant() const override { return 6.022140857E+23; }
  const double permitivity_vacuum() const override { return 8.854187817E-12; }
  const double elementary_charge() const override { return 1.6021766208E-19; }

  std::shared_ptr<PhysicalConstants> create(std::istream& istr) const override {
    return std::make_shared<CODATA2014>(istr); }
  void serialize(std::ostream& ostr) const override {
    ostr << class_name() << " "; }
  explicit CODATA2014(std::istream& istr);
  virtual ~CODATA2014() {}
};

inline std::shared_ptr<CODATA2014> MakeCODATA2014() {
  return std::make_shared<CODATA2014>();
}

/**
 * Physical constants from CODATA 2010
 * http://dx.doi.org/10.1103/RevModPhys.84.1527
 * http://dx.doi.org/10.1063/1.4724320
 */
class CODATA2010 : public PhysicalConstants {
 public:
  CODATA2010() : PhysicalConstants() {
    class_name_ = "CODATA2010";
    compute_derived_(); }

  const double boltzmann_constant() const override { return 1.3806488E-23; }
  const double avogadro_constant() const override { return 6.02214129E+23; }
  const double permitivity_vacuum() const override { return 8.854187817E-12; }
  const double elementary_charge() const override { return 1.602176565E-19; }

  std::shared_ptr<PhysicalConstants> create(std::istream& istr) const override {
    return std::make_shared<CODATA2010>(istr); }
  void serialize(std::ostream& ostr) const override {
    ostr << class_name() << " "; }
  explicit CODATA2010(std::istream& istr);
  virtual ~CODATA2010() {}
};

inline std::shared_ptr<CODATA2010> MakeCODATA2010() {
  return std::make_shared<CODATA2010>();
}

/**
 * Define custom physical constants.
 */
class PhysicalConstantsCustom : public PhysicalConstants {
 public:
  /**
    args:
    - boltzmann_constant
    - avogadro_constant
    - permitivity_vacuum
    - elementary_charge
   */
  explicit PhysicalConstantsCustom(argtype args = argtype());

  const double boltzmann_constant() const override {
    return boltzmann_constant_; }
  const double avogadro_constant() const override { return avogadro_constant_; }
  const double permitivity_vacuum() const override {
    return permitivity_vacuum_; }
  const double elementary_charge() const override { return elementary_charge_; }

  std::shared_ptr<PhysicalConstants> create(std::istream& istr) const override {
    return std::make_shared<PhysicalConstantsCustom>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit PhysicalConstantsCustom(std::istream& istr);
  virtual ~PhysicalConstantsCustom() {}

 private:
  double boltzmann_constant_;
  double avogadro_constant_;
  double permitivity_vacuum_;
  double elementary_charge_;
};

inline std::shared_ptr<PhysicalConstantsCustom> MakePhysicalConstantsCustom(
    argtype args = argtype()) {
  return std::make_shared<PhysicalConstantsCustom>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_PHYSICAL_CONSTANTS_H_
