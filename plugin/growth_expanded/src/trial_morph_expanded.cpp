#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "growth_expanded/include/compute_morph.h"
#include "growth_expanded/include/trial_morph.h"
#include "growth_expanded/include/trial_morph_expanded.h"

namespace feasst {

class MapTrialMorphExpanded {
 public:
  MapTrialMorphExpanded() {
    std::vector<std::vector<int> > seq = {{1}, {0}};
    auto obj = std::make_shared<TrialMorphExpanded>(seq);
    obj->deserialize_map()["TrialMorphExpanded"] = obj;
  }
};

static MapTrialMorphExpanded mapper_trial_morph_expanded_ = MapTrialMorphExpanded();

TrialMorphExpanded::TrialMorphExpanded(
    const std::vector<std::vector<int> > grow_seq,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialMorphExpanded";
  Arguments args_(args);
  args_.dont_check();
  bool add_previously = false;
  for (int step = 0;
       step < static_cast<int>(grow_seq.size());
       ++step) {
    const std::vector<int>& seq = grow_seq[step];
    argtype grow_args = args;
    argtype shrink_args = args;
    bool add = false;
    bool remove = false;
    int noskipgstage = 0;
    int noskipsstage = 0;
    for (int stage = 0; stage < static_cast<int>(seq.size()); ++stage) {
      // grow
      const int ptype0 = grow_seq[step][stage];
      if (ptype0 != -1) {
        int ptypem1 = -1;
        if (step > 0) ptypem1 = grow_seq[step-1][stage];
        if (ptypem1 == -1) {
          add = true;
          grow_args.insert({"particle_type" + str(noskipgstage), str(ptype0)});
        } else {
          ASSERT(!add, "cant have add and morph in same stage");
          grow_args.insert({"particle_type" + str(noskipgstage), str(ptypem1)});
          grow_args.insert({"particle_type_morph" + str(noskipgstage), str(ptype0)});
        }
        ++noskipgstage;
      }

      // shrink
      if (add_previously) {
        const int ptypem1 = grow_seq[step-1][stage];
        if (ptypem1 != -1) {
          remove = true;
          shrink_args.insert({"particle_type" + str(noskipsstage), str(ptypem1)});
          ++noskipsstage;
        }
      } else {
        int current_type, new_type;
        if (step == 0) {
          current_type = grow_seq[grow_seq.size() - 1][stage];
          new_type = grow_seq[grow_seq.size() - 2][stage];
        } else {
          ASSERT(step > 1, "step: " << step);
          current_type = grow_seq[step - 1][stage];
          new_type = grow_seq[step - 2][stage];
        }
        if (current_type != -1 && new_type != -1) {
          shrink_args.insert({"particle_type" + str(noskipsstage),
                              str(current_type)});
          shrink_args.insert({"particle_type_morph" + str(noskipsstage),
                              str(new_type)});
          ++noskipsstage;
        }
      }
    }
    if (add) {
      DEBUG(feasst_str(grow_args));
      grow_.push_back(MakeTrialAddMultiple(grow_args));
      DEBUG(feasst_str(grow_args));
      add = false;
      add_previously = true;
    } else {
      grow_.push_back(MakeTrialMorph(grow_args));
      add_previously = false;
    }
    if (remove) {
      DEBUG(feasst_str(shrink_args));
      shrink_.push_back(MakeTrialRemoveMultiple(shrink_args));
      remove = false;
    } else {
      shrink_.push_back(MakeTrialMorph(shrink_args));
    }
  }
  current_state_ = 0;
}

void TrialMorphExpanded::precompute(Criteria * criteria,
    System * system) {
  Trial::precompute(criteria, system);
  for (std::shared_ptr<Trial> gi : grow_) gi->precompute(criteria, system);
  for (std::shared_ptr<Trial> si : shrink_) si->precompute(criteria, system);
}

bool TrialMorphExpanded::attempt(Criteria * criteria,
    System * system,
    Random * random) {
  increment_num_attempts();
  const bool growing = random->coin_flip();
  DEBUG("growing: " << growing);
  bool accepted;
  if (growing) {
    accepted = grow_[current_state_]->attempt(criteria, system, random);
  } else {
    accepted = shrink_[current_state_]->attempt(criteria, system, random);
  }
  DEBUG("accepted: " << accepted);
  if (accepted) {
    increment_num_success_();
    if (growing) {
      ++current_state_;
    } else {
      --current_state_;
    }
    if (current_state_ == static_cast<int>(grow_.size())) {
      current_state_ = 0;
    } else if (current_state_ == -1) {
      current_state_ = static_cast<int>(grow_.size()) - 1;
    }
  }
  DEBUG("current_state: " << current_state_);
  return accepted;
}

std::shared_ptr<Trial> TrialMorphExpanded::create(std::istream& istr) const {
  return std::make_shared<TrialMorphExpanded>(istr);
}

TrialMorphExpanded::TrialMorphExpanded(std::istream& istr) : Trial(istr) {
  ASSERT(class_name_ == "TrialMorphExpanded", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2030 == version, "mismatch version: " << version);
  int dim1;
  istr >> dim1;
  grow_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      grow_[index] = grow_[index]->deserialize(istr);
    }
  }
  istr >> dim1;
  shrink_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      shrink_[index] = shrink_[index]->deserialize(istr);
    }
  }
  feasst_deserialize(&current_state_, istr);
}

void TrialMorphExpanded::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(2030, ostr);
  feasst_serialize_fstdr(grow_, ostr);
  feasst_serialize_fstdr(shrink_, ostr);
  feasst_serialize(current_state_, ostr);
}

}  // namespace feasst
