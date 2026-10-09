#include <algorithm>
#include <fstream>
#include "utils/include/utils.h" // find_in_list
#include "utils/include/file.h"  // remove_file
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "aniso/include/recursive_table.h"
#include "aniso/include/build_recursive_table.h"
#include "aniso/include/combine_parallel_files.h"

namespace feasst {

CombineParallelFiles::CombineParallelFiles(argtype * args) : Action(args) {
  class_name_ = "CombineParallelFiles";
  const std::string input_prefix = str("input_prefix", args);
  const std::string input_suffix = str("input_suffix", args, "");
  const std::string input2_suffix = str("input2_suffix", args, "");
  const std::string output_file = str("output_file", args, "");
  const int num_processors = integer("num_processors", args, -1);
  const std::string type = str("type", args);
  const bool recur_table = type.compare(0, 15, "recursive_table", 0, 15) == 0;
  RecursiveTable combined;
  if (recur_table) {
    for (int proc = 0; proc < num_processors; ++proc) {
      const std::string filename = input_prefix + str(proc) + input_suffix;
      DEBUG("filename " << filename);
      if (proc == 0) {
        DEBUG("proc0");
        combined = combined.from_file(filename);
        //combined.print_sizes();
      } else {
        RecursiveTable one;
        one = one.from_file(filename);
        //one.print_sizes();
        combined.add(one, proc, num_processors, type);
      }
    }
    DEBUG(combined.cutoff2d_.size());
  } else if (type == "contact_and_cutoff") {
    combined = combined.from_file(input_prefix + input_suffix);
    RecursiveTable cutoff;
    cutoff = cutoff.from_file(input_prefix + input2_suffix);
    combined.cutoff_ = cutoff.cutoff_;
    combined.cutoff2d_ = cutoff.cutoff2d_;
    combined.cutoff1d_ = cutoff.cutoff1d_;
  } else if (type == "contact_cutoff_energy") {
    combined = combined.from_file(input_prefix + input_suffix);
    RecursiveTable energy;
    energy = energy.from_file(input_prefix + input2_suffix);
    combined.energy_ = energy.energy_;
    combined.energy3d_ = energy.energy3d_;
    combined.energy2d_ = energy.energy2d_;
  } else if (type == "trained") {
    std::vector<std::vector<rbin> > rbins_in_proc(num_processors);
    rbin max_bins;
    for (int proc = 0; proc < num_processors; ++proc) {
      const std::string filename = input_prefix + str(proc) + input_suffix;
      rbins_in_proc[proc] = read_rbins(filename);
    }

    // union rbins with the same elements but possibly different criteria
    std::vector<rbin> rbins;
    for (int proc = 0; proc < num_processors; ++proc) {
      for (const rbin& rb : rbins_in_proc[proc]) {
        int index;
        if (find_in_list(rb.first, rbins, &index)) {
          if (rbins[index].second < rb.second) {
            rbins[index].second = rb.second;
          }
        } else {
          rbins.push_back(rb);
        }
      }
    }
    sort_and_print(&rbins, input_prefix+input_suffix);
  } else {
    FATAL("unrecognized file type: " << type);
  }

  if (recur_table || type == "contact_and_cutoff" || type == "contact_cutoff_energy") {
    std::stringstream ss;
    combined.serialize(ss);
    std::string filename;
    if (!output_file.empty()) {
      filename = output_file;
    } else {
      if (recur_table) {
        filename = input_prefix+input_suffix;
      } else {
        filename = input_prefix+input_suffix+input2_suffix;
      }
    }
    std::ofstream file(filename);
    ASSERT(file.good(), "Error handing filename: " << filename);
    file << ss.str();
  }

  // If successful, remove input files
  if (boolean("cleanup", args, true)) {
    if (recur_table || type == "trained") {
      DEBUG("cleaning");
      for (int proc = 0; proc < num_processors; ++proc) {
        const std::string filename = input_prefix + str(proc) + input_suffix;
        remove_file(filename);
      }
    } else if (type == "contact_and_cutoff" || type == "contact_cutoff_energy") {
      remove_file(input_prefix + input_suffix);
      remove_file(input_prefix + input2_suffix);
    }
  }
}
CombineParallelFiles::CombineParallelFiles(argtype args) : Action(&args) {
  feasst_check_all_used(args);
}
CombineParallelFiles::~CombineParallelFiles() {}

FEASST_MAPPER(CombineParallelFiles,);

CombineParallelFiles::CombineParallelFiles(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 4367 && version <= 4367, "mismatch version: " << version);
}

void CombineParallelFiles::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(4367, ostr);
}

}  // namespace feasst
