#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/backmap.h"

namespace feasst {

FEASST_MAPPER(Backmap, argtype({{"output_file", "place_holder"}}));

Backmap::Backmap(argtype * args) : AnalyzeWriteOnly(args) {
  set_append();
  ASSERT(!output_file().empty(), "file name is required");

  if (used("site0", *args)) {
    WARN("site[i] is a deprecated argument. Please use sites instead.");
    DEBUG("parse site types to backmap");
    int type = 0;
    std::stringstream key;
    key << "site" << type;
    while (used(key.str(), *args)) {
      site_types_.push_back(integer(key.str(), args));
      key.str("");
      key << "fstprt" << type;
      site_fstprt_.push_back(str(key.str(), args));
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << "site" << type;
    }
  } else {
    site_type_names_ = split(feasst::str("sites", args, ""), ',');
    site_fstprt_ = split(feasst::str("fstprts", args, ""), ',');
    ASSERT(site_type_names_.size() == site_fstprt_.size(), "Number of sites:" <<
      site_type_names_.size() << " != number of fstprts:" << site_fstprt_.size());
    DEBUG(feasst_str(site_type_names_));
  }
  DEBUG(feasst_str(site_types_));
  DEBUG(feasst_str(site_fstprt_));

  args->insert({"append", "true"}); // always append
  xyz_ = FileXYZ(args);
  vmd_ = FileVMD(args);
}
Backmap::Backmap(argtype args) : Backmap(&args) {
  feasst_check_all_used(args);
}

void Backmap::add_backmap_particles_() {
  argtype config_args = {{"wrap", "false"}};
  config_args.insert({"particle_type",feasst_str(site_fstprt_)});
//  int type = 0;
//  for (const std::string& fstprt : site_fstprt_) {
//    config_args.insert({"particle_type"+str(type), fstprt});
//    ++type;
//  }
  all_atom_ = MakeConfiguration(config_args);
}

void Backmap::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  const std::string name = output_file(mc->criteria());
  ASSERT(!name.empty(), "file name required. Did you forget to " <<
    "Analyze::set_output_file()?");
  const Configuration& config = configuration(mc->system());
  if (site_type_names_.size() > 0) {
    for (const std::string& name : site_type_names_) {
      site_types_.push_back(config.site_type_name_to_index(name));
    }
  }
  for (const int site_type : site_types_) {
    const int part_type = config.site_type_to_particle_type(site_type);
    ASSERT(config.unique_type(part_type, site_type).is_anisotropic(),
      "site_type:" << site_type << " is not anisotropic.");
  }

  add_backmap_particles_();

  // write vmd
  std::stringstream ss;
  ss << name << ".vmd";
  vmd_.write(ss.str(), *all_atom_, name);

  // write immediately like Movie does so trajectories are 1-1
  write(*mc);
}

std::string Backmap::write(const MonteCarlo& mc) {
  const Configuration& orig_config = configuration(mc.system());
  add_backmap_particles_();
  all_atom_->set(std::make_shared<Domain>(orig_config.domain()));
  const Select& all = orig_config.selection_of_all();
  std::vector<std::shared_ptr<TrialSelectParticle> > sels;
  for (int site_type : site_types_) {
    sels.push_back(
      MakeTrialSelectParticle({{"particle_type", str(site_type)}}));
  }
  PerturbRotate rotate;
  PerturbTranslate translate;
  Position origin;
  const int dimen = orig_config.dimension();
  origin.set_to_origin(dimen);
  RotationMatrix rot_mat;
  rot_mat.set_size(dimen);
  System sys;
  sys.add(all_atom_);
  Configuration * new_config = sys.get_configuration();
  for (int select_index = 0;
       select_index < all.num_particles();
       ++select_index) {
    const int part_index = all.particle_index(select_index);
    const Particle& part = orig_config.select_particle(part_index);
    //const int part_type = part.type();
    for (int site_index : all.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      const int site_type = site.type();
      int findex;
      if (find_in_list(site_type, site_types_, &findex)) {
        std::shared_ptr<TrialSelectParticle> sel = sels[findex];
        const Position& center = site.position();
        const Euler& euler = site.euler();
        new_config->add_particle_of_type(findex);
        DEBUG("euler: " << euler.str());
        euler.compute_rotation_matrix(&rot_mat);
        DEBUG("rot_mat " << rot_mat.str());
        sel->select_particle(new_config->num_particles()-1, *new_config);
        DEBUG("size pos:" << site.position().str());
        DEBUG("sel: " << sel->mobile().str());
        DEBUG("sel prop size:" << sel->mobile().site_properties().size());
        DEBUG("num particles:" << new_config->num_particles());
        sel->set_mobile_original(&sys);
        rotate.move(origin, rot_mat, &sys, sel.get());
        translate.move(center, &sys, sel.get());
      }
    }
  }
  xyz_.write(output_file(mc.criteria()), *new_config);
  return std::string("");
}

void Backmap::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2634, ostr);
  feasst_serialize(site_types_, ostr);
  feasst_serialize_fstobj(xyz_, ostr);
  feasst_serialize_fstobj(vmd_, ostr);
}

Backmap::Backmap(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2634, "version mismatch:" << version);
  feasst_deserialize(&site_types_, istr);
  feasst_deserialize_fstobj(&xyz_, istr);
  feasst_deserialize_fstobj(&vmd_, istr);
}

}  // namespace feasst
