#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "morph/include/perturb_particle_type.h"
#include "morph/include/compute_morph.h"
#include "morph/include/trial_morph.h"

namespace feasst {

FEASST_MAPPER(TrialMorph, argtype({{"particle_type", "0"},
                                   {"particle_type_morph", "1"}}));

TrialMorph::TrialMorph(argtype * args) : Trial(args) {
  class_name_ = "TrialMorph";
  set_description("TrialMorph");
  argtype stage_args = get_stage_args(args);
  if (used("particle_type", *args)) {
    std::vector<std::string> pts = split(str("particle_type", args), ',');
    std::vector<std::string> ptms = split(str("particle_type_morph", args), ',');
    ASSERT(pts.size() == ptms.size(), "The number of particle_type:" <<
      pts.size() << " should match the number of particle_type_more:" <<
      ptms.size());
    for (int ipt = 0; ipt < static_cast<int>(pts.size()); ++ipt) {
      const std::string& pt = pts[ipt];
      const std::string& ptm = ptms[ipt];
      ASSERT(pt != ptm,
        "Attempting to morph particles into the same type:" << pt);
      auto sel = MakeTrialSelectParticle({{"particle_type", pt},
                                          {"exclude_perturbed", "true"}});
      argtype tmp_stage_args = stage_args;
      add_stage(sel, MakePerturbParticleType({{"type", ptm}}), &tmp_stage_args);
    }
  } else {
    std::string start("particle_type");
    std::stringstream key;
    int type = 0;
    key << start << type;
    //INFO("args: " << str(*args));
    while (used(key.str(), *args)) {
      const std::string type_str = feasst::str(type);
      const std::string ptype = str(key.str(), args);
      const std::string pmt = str("particle_type_morph" + type_str, args);
      ASSERT(ptype != pmt, "Attempting to morph particles into the same type: " <<
        ptype);
      auto sel = MakeTrialSelectParticle({{"particle_type", ptype},
                                          {"exclude_perturbed", "true"}});
      argtype tmp_stage_args = stage_args;
      add_stage(sel, MakePerturbParticleType({{"type", pmt}}), &tmp_stage_args);
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
    ASSERT(type > 0, "required arguments not used");
    WARN("Deprecated TrialMorph::particle_type[i]->particle_type");
  }
  set(MakeComputeMorph());
}
TrialMorph::TrialMorph(argtype args) : TrialMorph(&args) {
  feasst_check_all_used(args);
}

TrialMorph::TrialMorph(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8967, "mismatch version: " << version);
}

void TrialMorph::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(8967, ostr);
}

}  // namespace feasst
