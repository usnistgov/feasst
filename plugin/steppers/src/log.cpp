#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/log.h"

namespace feasst {

FEASST_MAPPER(Log,);

Log::Log(argtype args) : Log(&args) { feasst_check_all_used(args); }
Log::Log(argtype * args) : AnalyzeWriteOnly(args) {
  if (boolean("append", args, true)) {
    set_append();
  } else {
    ERROR("append is required");
  }
  max_precision_ = boolean("max_precision", args, false);
  include_bonds_ = boolean("include_bonds", args, true);
  format_ = str("format", args, "csv");
  ASSERT(format_ == "csv" || format_ == "vertical",
    "unrecognized format: " << format_);
}

void Log::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string Log::header(const MonteCarlo& mc) const {
  const System& system = mc.system();
  std::stringstream ss;
  if (format_ == "vertical") {
    ss << "#";
  }
  ss << system.status_header();
  ss << ",";
  if (mc.is_criteria_set()) {
    ss << mc.criteria().status_header(system, include_bonds_);
    ss << ",";
    // print number of trials here instead of TrialFactory header because
    // multiple factories makes it redundant.
    ss << ",trial"
       << mc.trial_factory().status_header()
       << std::endl;
  }
  return ss.str();
}

std::string Log::write(const MonteCarlo& mc) {
  const System& system = mc.system();
  const TrialFactory& trial_factory = mc.trial_factory();
  // ensure the following order matches the header from initialization.
  std::stringstream ss;
//  if (format_ == "csv") {
//    if (rewrite_header()) {
//      ss << header(mc);
//    }
//  }
  std::stringstream values;
  values << system.status();
  values << ",";
  if (mc.is_criteria_set()) {
    values << mc.criteria().status(system, max_precision_, include_bonds_, &bond_visitor_);
    values << ",";
    values << "," << trial_factory.num_attempts()
       << trial_factory.status()
       << std::endl;
  }
  if (format_ == "csv") {
    ss << values.str();
  } else if (format_ == "vertical") {
    std::string hdr = header(mc);
    // # remove # character added to header, then remove end lines
    hdr.erase(0,1);
    auto hline = split(hdr, ',');
    if (hline.back().back() == '\n') {
      hline.back().pop_back();
    }
    auto vline = split(values.str(), ',');
    if (vline.back().back() == '\n') {
      vline.back().pop_back();
    }
    const int hsize = static_cast<int>(hline.size());
    const int vsize = static_cast<int>(vline.size());
    if (hsize != vsize) {
      WARN("header has " << hsize << " elements, but there are " << vsize <<
        " values.");
    }
    for (int index = 0; index < std::max(hsize, vsize); ++index) {
      if (index >= hsize) {
        ss << "NaN";
      } else {
        ss << hline[index];
      }
      ss << "=";
      if (index >= vsize) {
        ss << "NaN";
      } else {
        ss << vline[index];
      }
      ss << std::endl;
    }
    ss << std::endl;
  }
  return ss.str();
}

void Log::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(669, ostr);
  feasst_serialize(max_precision_, ostr);
  feasst_serialize(include_bonds_, ostr);
  feasst_serialize_fstobj(bond_visitor_, ostr);
  feasst_serialize(format_, ostr);
}

Log::Log(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 668 && version <= 669, "version mismatch:" << version);
  feasst_deserialize(&max_precision_, istr);
  feasst_deserialize(&include_bonds_, istr);
  feasst_deserialize_fstobj(&bond_visitor_, istr);
  if (version >= 669) {
    feasst_deserialize(&format_, istr);
  }
}

}  // namespace feasst
