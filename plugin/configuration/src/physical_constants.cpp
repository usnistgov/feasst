#include <cmath>
#include "utils/include/serialize_extra.h"
#include "utils/include/arguments.h"
#include "math/include/constants.h"
#include "configuration/include/physical_constants.h"

namespace feasst {

std::map<std::string, std::shared_ptr<PhysicalConstants> >&
    PhysicalConstants::deserialize_map() {
  static std::map<std::string, std::shared_ptr<PhysicalConstants> >* ans =
     new std::map<std::string, std::shared_ptr<PhysicalConstants> >();
  return *ans;
}

std::shared_ptr<PhysicalConstants> PhysicalConstants::deserialize(
    std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void PhysicalConstants::compute_derived_() {
  ideal_gas_constant_ = boltzmann_constant()*avogadro_constant();
  charge_conversion_ = std::pow(elementary_charge(), 2)/
  (4*PI*permitivity_vacuum()*1e3/1e10/avogadro_constant());
}

FEASST_MAPPER(CODATA2018,);

CODATA2018::CODATA2018(std::istream& istr) : PhysicalConstants(istr) {
  compute_derived_();
}

FEASST_MAPPER(CODATA2014,);

CODATA2014::CODATA2014(std::istream& istr) : PhysicalConstants(istr) {
  compute_derived_();
}

FEASST_MAPPER(CODATA2010,);

CODATA2010::CODATA2010(std::istream& istr) : PhysicalConstants(istr) {
  compute_derived_();
}

FEASST_MAPPER(PhysicalConstantsCustom, argtype(
  {{"boltzmann_constant", "0"},
   {"avogadro_constant", "0"},
   {"permitivity_vacuum", "0"},
   {"elementary_charge", "0"}}));

PhysicalConstantsCustom::PhysicalConstantsCustom(argtype args)
  : PhysicalConstants() {
  class_name_ = "PhysicalConstantsCustom";
  boltzmann_constant_ = dble("boltzmann_constant", &args);
  avogadro_constant_ = dble("avogadro_constant", &args);
  permitivity_vacuum_ = dble("permitivity_vacuum", &args);
  elementary_charge_ = dble("elementary_charge", &args);
  compute_derived_();
  feasst_check_all_used(args);
}

void PhysicalConstantsCustom::serialize(std::ostream& ostr) const {
  ostr << class_name() << " ";
  feasst_serialize(boltzmann_constant_, ostr);
  feasst_serialize(avogadro_constant_, ostr);
  feasst_serialize(permitivity_vacuum_, ostr);
  feasst_serialize(elementary_charge_, ostr);
}

PhysicalConstantsCustom::PhysicalConstantsCustom(std::istream& istr)
  : PhysicalConstants(istr) {
  feasst_deserialize(&boltzmann_constant_, istr);
  feasst_deserialize(&avogadro_constant_, istr);
  feasst_deserialize(&permitivity_vacuum_, istr);
  feasst_deserialize(&elementary_charge_, istr);
  compute_derived_();
}

void CODATA2018::serialize(std::ostream& ostr) const {
  ostr << class_name() << " "; }
void CODATA2014::serialize(std::ostream& ostr) const {
  ostr << class_name() << " "; }
void CODATA2010::serialize(std::ostream& ostr) const {
  ostr << class_name() << " "; }

}  // namespace feasst
