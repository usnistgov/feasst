#include <vector>
#include <string>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "system/include/model.h"
#include "system/include/visit_model.h"
#include "system/include/potential.h"
#include "system/include/bond_visitor.h"
#include "system/include/potential_factory.h"

namespace feasst {

PotentialFactory::PotentialFactory() {
  initialize_bond_visitor_();
}
PotentialFactory::~PotentialFactory() {}

void PotentialFactory::initialize_bond_visitor_() {
  if (static_cast<int>(potentials_.size()) == 1) {
    if (potentials_[0]->visit_model().class_name() == "DontVisitModel") {
      bonds_.clear();
    }
  } else {
    if (static_cast<int>(bonds_.size()) == 0) {
      bonds_.push_back(std::make_shared<BondVisitor>());
    }
  }
}

void PotentialFactory::add(std::shared_ptr<Potential> potential) {
  potentials_.push_back(potential);
  initialize_bond_visitor_();
}

void PotentialFactory::set(const int index,
                           std::shared_ptr<Potential> potential) {
  potentials_[index] = potential;
  initialize_bond_visitor_();
}

void PotentialFactory::precompute(Configuration * config) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    DEBUG("precomputing " << potential->model().class_name());
    potential->precompute(config);
  }
}

void PotentialFactory::precompute(const int index, Configuration * config) {
  DEBUG("precomputing index: " << index);
  potentials_[index]->precompute(config);
}

double PotentialFactory::energy(Configuration * config) {
  double en = 0;
  int index = 0;
  while ((index < num()) && (opt_overlap_ == 0 || (en < NEAR_INFINITY/10.))) {
    const double potential_en = potentials_[index]->energy(config);
    DEBUG("potential index: " << index << " potential energy: " << potential_en);
    en += potential_en;
    ++index;
  }
  DEBUG("en " << en);
  DEBUG(str());
  double bond_en = 0;
  for (std::shared_ptr<BondVisitor> bn : bonds_) {
    bn->compute_all(*config);
    bond_en += bn->energy();
  }
  DEBUG("bond_en " << bond_en);
  return en + bond_en;
}

double PotentialFactory::select_energy(const Select& select, Configuration * config) {
  double en = 0;
  int index = 0;
  //while (index < static_cast<int>(potentials_.size())) {
  while ((index < static_cast<int>(potentials_.size())) &&
         (opt_overlap_ == 0 || (en < NEAR_INFINITY))) {
    DEBUG("index " << index);
    en += potentials_[index]->select_energy(select, config);
    ++index;
  }
  DEBUG("en " << en);
  DEBUG(str());
  ASSERT(!std::isinf(en), "en: " << en << " is inf.");
  ASSERT(!std::isnan(en), "en: " << en << " is nan.");
  double bond_en = 0;
  for (std::shared_ptr<BondVisitor> bn : bonds_) {
    bn->compute_all(select, *config);
    bond_en += bn->energy();
  }
  DEBUG("bond en " << bond_en);
  ASSERT(!std::isinf(bond_en), "bond_en: " << bond_en << " is inf.");
  ASSERT(!std::isnan(bond_en), "bond_en: " << bond_en << " is nan.");
  return en + bond_en;
}

std::vector<double> PotentialFactory::stored_energy_profile() const {
  std::vector<double> en;
  for (const std::shared_ptr<Potential>& potential : potentials_) {
    en.push_back(potential->stored_energy());
  }
  return en;
}

double PotentialFactory::stored_energy() const {
  std::vector<double> en = stored_energy_profile();
  return std::accumulate(en.begin(), en.end(), 0.);
}

std::string PotentialFactory::str() const {
  std::stringstream ss;
  ss << "PotentialFactory: ";
  for (const std::shared_ptr<Potential>& potential : potentials_) {
    ss << potential->stored_energy() << " ";
  }
  return ss.str();
}

void PotentialFactory::revert(const Select& select) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    potential->revert(select);
  }
}

void PotentialFactory::finalize(const Select& select, Configuration * config) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    potential->finalize(select, config);
  }
}

void PotentialFactory::serialize(std::ostream& sstr) const {
  feasst_serialize_version(8656, sstr);
  feasst_serialize(potentials_, sstr);
  feasst_serialize(opt_overlap_, sstr);
  feasst_serialize(user_name_, sstr);
  feasst_serialize(bonds_, sstr);
}

PotentialFactory::PotentialFactory(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version >= 8655 && version <= 8656, "unrecognized verison: " << version);
  // HWH for unknown reasons, this does not work
  // feasst_deserialize(&potentials_, sstr);
  int dim1;
  sstr >> dim1;
  potentials_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    sstr >> existing;
    if (existing != 0) {
      potentials_[index] = std::make_shared<Potential>(sstr);
    }
  }
  feasst_deserialize(&opt_overlap_, sstr);
  if (version >= 8656) {
    feasst_deserialize(&user_name_, sstr);
    //  HWH for unknown reasons, this function template does not work.
    //feasst_deserialize(&bonds_, sstr);
    {
      int dim1;
      sstr >> dim1;
      bonds_.resize(dim1);
      for (int index = 0; index < dim1; ++index) {
        //feasst_deserialize((*vector)[index], istr);
        int existing;
        sstr >> existing;
        if (existing != 0) {
          bonds_[index] = std::make_shared<BondVisitor>(sstr);
        }
      }
    }
  }
}

void PotentialFactory::load_cache(const bool load) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    potential->load_cache(load);
  }
}

void PotentialFactory::unload_cache(const PotentialFactory& factory) {
  ASSERT(num() == factory.num(), "size mismatch");
  for (int ip = 0; ip < num(); ++ip) {
    potentials_[ip]->unload_cache(*factory.potentials_[ip]);
  }
}

void PotentialFactory::check(const Configuration& config) const {
  for (const std::shared_ptr<Potential>& potential : potentials_) {
    potential->check(config);
  }
}

void PotentialFactory::synchronize_(const PotentialFactory& factory,
  const Select& perturbed) {
  for (int index = 0; index < num(); ++index) {
    potentials_[index]->synchronize_(*factory.potentials()[index], perturbed);
  }
}

void PotentialFactory::change_volume(const double delta_volume,
    const int dimension, Configuration * config) {
  for (int index = 0; index < num(); ++index) {
    potentials_[index]->change_volume(delta_volume, dimension, config);
  }
}

void PotentialFactory::remove_opt_overlap() {
  opt_overlap_ = 0;
}

}  // namespace feasst
