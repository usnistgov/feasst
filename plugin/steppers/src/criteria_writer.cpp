#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/criteria_writer.h"

namespace feasst {

CriteriaWriter::CriteriaWriter(argtype * args) : AnalyzeWriteOnly(args) {}
CriteriaWriter::CriteriaWriter(argtype args) : CriteriaWriter(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(CriteriaWriter,);

std::string CriteriaWriter::write(const MonteCarlo& mc) {
  // ensure the following order matches the header from initialization.
  std::stringstream ss;
  ss << mc.criteria().write() << std::endl;
  return ss.str();
}

void CriteriaWriter::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(993, ostr);
}

CriteriaWriter::CriteriaWriter(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 993, "version mismatch:" << version);
}

}  // namespace feasst
