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
  const int num_processors = integer("num_processors", args);
  const std::string type = str("type", args);
  if (type == "recursive_table") {
    RecursiveTable reader;
    std::vector<RecursiveTable> tables;
    //const std::vector<int> bins = {0,0,0,0,0};  // temporary for testing
    for (int proc = 0; proc < num_processors; ++proc) {
      const std::string filename = input_prefix + str(proc) + input_suffix;
      DEBUG("filename " << filename);
      tables.push_back(reader.from_file(filename));
      DEBUG(proc << " " << tables[proc].cutoff2d_.size());
    }
    RecursiveTable combined = reader.from_parallel_build(tables);
    DEBUG(combined.cutoff2d_.size());
//    INFO(tables[0].contact_[0][0].nested(bins).data()[0][0][0][0][0]);
    //INFO( combined.contact_[0][0].nested(bins).data()[0][0][0][0][1]);
    std::stringstream ss;
    combined.serialize(ss);
    const std::string filename = input_prefix+input_suffix;
    std::ofstream file(filename);
    ASSERT(file.good(), "Error handing filename: " << filename);
    file << ss.str();
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

  // If successful, remove input files
  if (boolean("cleanup", args, true)) {
    DEBUG("cleaning");
    for (int proc = 0; proc < num_processors; ++proc) {
      const std::string filename = input_prefix + str(proc) + input_suffix;
      remove_file(filename);
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
