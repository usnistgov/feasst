#include <gtest/gtest.h>
#include "chain/include/perturb_regrow.h"
#include "system/include/visit_model_intra.h"
#include "system/include/model_hard_sphere.h"
#include "configuration/include/file_xyz.h"

namespace feasst {

//TEST(PerturbRegrow, regrow) {
//  seed_random_by_date();
//  System system;
//  { Configuration config;
//    config.set_domain(Domain().set_cubic(8));
//    config.add_particle_type("../forcefield/data.chain10");
//    config.add_particle_of_type(0);
//    system.add(config); }
//
//  { Potential potential;
//    potential.set_model(std::make_shared<ModelHardSphere>());
//    auto visitor = std::make_shared<VisitModelIntra>();
//    visitor->set_intra_cut(1);
//    potential.set_visit_model(visitor);
//    system.add_to_unoptimized(potential); }
//
//  PerturbRegrow regrow;
//  regrow.select(&system);
//  //regrow.select_random_end_segment_in_particle(0, system.configuration());
//  regrow.perturb(&system);
//  FileXYZ().write_for_vmd("tmp/regrow.xyz", system.configuration());
//}

}  // namespace feasst
