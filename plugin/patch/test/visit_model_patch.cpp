#include <gtest/gtest.h>
#include "patch/include/visit_model_patch.h"
#include "core/include/file_xyz.h"
#include "core/include/model_square_well.h"
#include "core/include/perturb_translate.h"
#include "core/include/perturb_rotate.h"
#include "core/include/file_xyz.h"

namespace feasst {

TEST(VisitModelPatch, patch_one) {
  Configuration config;
  config.set_domain(Domain().set_cubic(10.));
  config.add_particle_type("../plugin/patch/forcefield/data.patch_one");
  config.set_model_param("cutoff", 0, 3.);
  config.set_model_param("cutoff", 1, 3.);
  FileXYZ().load("../plugin/patch/test/data/patch5.xyz", &config);
  ModelSquareWell model;
  VisitModelPatch visitor;
  config.add(Group().add_site_type(0));
  visitor.compute(model, &config, 1);
  EXPECT_NEAR(-3., visitor.energy(), NEAR_ZERO);
}

TEST(VisitModelPatch, patch_one_2body) {
  System system;
  { Configuration config;
    config.set_domain(Domain().set_cubic(10.));
    auto angle = std::make_shared<ModelParam>();
    angle->set_name("angle");
    config.add(angle);
    config.add_particle_type("../plugin/patch/forcefield/data.patch_one");
    EXPECT_EQ(0, config.model_params().select("angle")->value(0));
    EXPECT_EQ(5, config.model_params().select("angle")->value(1));
    config.add(Group().add_site_type(0));
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    system.add(config); }

  { Potential potential;
    potential.set_model(std::make_shared<ModelSquareWell>());
    auto visitor = std::make_shared<VisitModelPatch>();
    visitor->cpa_sq_ = pow(cos(5./180.*PI), 2);
    potential.set_visit_model(visitor);
    potential.set_group_index(1); // optimization: loop through centers only.
    system.add(potential); }

  PerturbTranslate trans;
  trans.set_selection(SelectList().particle(1, system.configuration()));
  auto traj = Position().set_vector({1.4, 0., 0.});
  trans.translate_selection(traj, &system);

  PerturbRotate rotate;
  rotate.set_selection(SelectList().particle(1, system.configuration()));
  auto axis = Position().set_vector({0., 0., 1.});
  FileXYZ file;
  file.write_for_vmd("tmp/patch.xyz", system.configuration());
  rotate.rotate_selection(traj,
    RotationMatrix().axis_angle(axis, 180.),
    &system);
  file.set_append(1);
  file.write("tmp/patch.xyz", system.configuration());

  EXPECT_NEAR(-1., system.energy(), NEAR_ZERO);

  rotate.revert();
  rotate.rotate_selection(traj,
    RotationMatrix().axis_angle(axis, 180 - 4.999),
    &system);
  EXPECT_NEAR(-1., system.energy(), NEAR_ZERO);

  rotate.revert();
  rotate.rotate_selection(traj,
    RotationMatrix().axis_angle(axis, 180 - 5.001),
    &system);
  EXPECT_NEAR(0., system.energy(), NEAR_ZERO);
}

} // namespace feasst

