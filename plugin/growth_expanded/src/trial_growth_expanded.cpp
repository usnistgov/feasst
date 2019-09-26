#include "growth_expanded/include/trial_growth_expanded.h"
#include "monte_carlo/include/trial_add.h"

namespace feasst {

class MapTrialComputeGrowAdd {
 public:
  MapTrialComputeGrowAdd() {
    auto obj = std::make_shared<TrialComputeGrowAdd>();
    obj->deserialize_map()["TrialComputeGrowAdd"] = obj;
  }
};

static MapTrialComputeGrowAdd mapper_trial_compute_grow_add_ = MapTrialComputeGrowAdd();

std::shared_ptr<TrialCompute> TrialComputeGrowAdd::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGrowAdd>(istr);
}

TrialComputeGrowAdd::TrialComputeGrowAdd(std::istream& istr)
  : TrialCompute(istr) {
  ASSERT(class_name_ == "TrialComputeGrowAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(154 == version, "mismatch version: " << version);
}


void TrialComputeGrowAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(154, ostr);
}

class MapTrialComputeGrowRemove {
 public:
  MapTrialComputeGrowRemove() {
    auto obj = std::make_shared<TrialComputeGrowRemove>();
    obj->deserialize_map()["TrialComputeGrowRemove"] = obj;
  }
};

static MapTrialComputeGrowRemove mapper_trial_compute_grow_remove_ = MapTrialComputeGrowRemove();

std::shared_ptr<TrialCompute> TrialComputeGrowRemove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGrowRemove>(istr);
}

TrialComputeGrowRemove::TrialComputeGrowRemove(std::istream& istr)
  : TrialCompute(istr) {
  ASSERT(class_name_ == "TrialComputeGrowRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(225 == version, "mismatch version: " << version);
}


void TrialComputeGrowRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(225, ostr);
}

class MapTrialComputeGrow {
 public:
  MapTrialComputeGrow() {
    auto obj = std::make_shared<TrialComputeGrow>();
    obj->deserialize_map()["TrialComputeGrow"] = obj;
  }
};

static MapTrialComputeGrow mapper_trial_compute_grow_ = MapTrialComputeGrow();

std::shared_ptr<TrialCompute> TrialComputeGrow::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGrow>(istr);
}

TrialComputeGrow::TrialComputeGrow(std::istream& istr)
  : TrialCompute(istr) {
  ASSERT(class_name_ == "TrialComputeGrow", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(155 == version, "mismatch version: " << version);
}


void TrialComputeGrow::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(155, ostr);
}

class MapTrialGrowthExpanded {
 public:
  MapTrialGrowthExpanded() {
    auto obj = std::make_shared<TrialGrowthExpanded>(MakeTrialAdd(
      {{"particle_type", "0"}}), MakeTrialAdd({{"particle_type", "0"}}));
    obj->deserialize_map()["TrialGrowthExpanded"] = obj;
  }
};

static MapTrialGrowthExpanded mapper_trial_growth_expanded_ = MapTrialGrowthExpanded();

std::shared_ptr<Trial> TrialGrowthExpanded::create(std::istream& istr) const {
  return std::make_shared<TrialGrowthExpanded>(istr);
}

TrialGrowthExpanded::TrialGrowthExpanded(std::istream& istr)
  : Trial(istr) {
  ASSERT(class_name_ == "TrialGrowthExpanded", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(793 == version, "mismatch version: " << version);
}


void TrialGrowthExpanded::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(793, ostr);
}

}  // namespace feasst
