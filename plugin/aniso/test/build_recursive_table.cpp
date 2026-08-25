#include "utils/test/utils.h"
#include "utils/include/file.h"  // compare_files
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/recursive_table.h"
#include "aniso/include/build_recursive_table.h"

namespace feasst {

TEST(BuildRecursiveTable, serialize) {
  BuildRecursiveTable obj;
  BuildRecursiveTable obj2 = test_serialize(obj);
}

TEST(BuildRecursiveTable, parallel_VERY_LONG) {
  const int num_proc = 3;
  const argtype config_args = {{"side_length", "50,50,50"}, {"periodic", "false,false,false"}, {"particle_type", "t1:../plugin/aniso/particle/aniso_l10_patch.txt"}, {"add_num_t1_particles", "2"}};

  // mayer sampling training data
  argtype mayer_config_args = config_args;
  mayer_config_args.insert({"group", "fixed,mobile"});
  mayer_config_args.insert({"fixed_particle_index", "0"});
  mayer_config_args.insert({"mobile_particle_index", "1"});
  mayer_config_args.insert({"set_cutoff_min_to_sigma", "true"});
  const std::string contact_only = "true";

  if (false) {
  //if (true) {
    WARN("comment above to generate necessary files");
  } else {

  MakeMonteCarlo({{
    {"Configuration", mayer_config_args},
    {"Potential", {{"Model", "HardSphere"}}},
    {"RefPotential", {{"ref", "hs"}, {"Model", "HardSphere"}, {"sigma", "0"}, {"sigmaA_A", "1"}}},
    {"ThermoParams", {{"beta", "1"}}},
    {"MayerSampling", {{"trials_per_cycle", "1e4"}, {"cycles_to_complete", "1e2"}}},
    {"TrialTranslate", {{"new_only", "true"}, {"ref", "hs"}, {"tunable_param", "1"}, {"group", "mobile"}}},
    {"TrialRotate", {{"new_only", "true"}, {"ref", "hs"}, {"tunable_param", "1"}, {"group", "mobile"}}},
    {"Tune", {{}}},
    {"Run", {{"until", "complete"}}},
    {"Remove", {{"name", "Tune"}}}, //,CriteriaWriter,Log,Movie
    {"CriteriaWriter", {{"trials_per_write", "1e4"}, {"output_file", "tmp/recurb2.csv"}}},
    {"MayerSampling", {{"trials_per_cycle", "1e4"}, {"cycles_to_complete", "1e2"}, {"training_file", "tmp/recur_training.csv"}}},
    {"Run", {{"until", "complete"}}},
  }}, true);
}

  // build base
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "HardSphere"}}},
    {"BuildRecursiveTable", {{"stage","base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/base_serial"}, {"contact_only", contact_only}}},
  }}, true);
  RecursiveTable base_serial = RecursiveTable().from_file("tmp/base_serial");
  for (int proc = 0; proc < num_proc; ++proc) {
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "HardSphere"}}},
      {"BuildRecursiveTable", {{"stage", "base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/base"+str(proc)}, {"contact_only", contact_only}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/base"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  RecursiveTable base_parallel = RecursiveTable().from_file("tmp/base");
  EXPECT_TRUE(base_serial.is_equal(base_parallel));

  // train
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "HardSphere"}}},
    {"BuildRecursiveTable", {{"stage", "train"}, {"mayer_training_file", "tmp/recur_training.csv"}, {"base_table_file", "tmp/base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/train_serial"}, {"contact_only", contact_only}, {"min_criteria", "5"}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "HardSphere"}}},
      {"BuildRecursiveTable", {{"stage", "train"}, {"mayer_training_file", "tmp/recur_training.csv"}, {"base_table_file", "tmp/base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/train"+str(proc)}, {"contact_only", contact_only}, {"min_criteria", "5"}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "trained"}, {"input_prefix", "tmp/train"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  EXPECT_TRUE(compare_files("tmp/train", "tmp/train_serial"));

  // build nested tables using the base table and train file
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "HardSphere"}}},
    {"BuildRecursiveTable", {{"stage", "build_contact"}, {"base_table_file", "tmp/base"}, {"trained_file", "tmp/train"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/recur_serial"}, {"contact_only", contact_only}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    INFO("proc " << proc);
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "HardSphere"}}},
      {"BuildRecursiveTable", {{"stage", "build_contact"}, {"base_table_file", "tmp/base"}, {"trained_file", "tmp/train"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/recur"+str(proc)}, {"contact_only", contact_only}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
    INFO("combining");
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/recur"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  RecursiveTable recur_serial = RecursiveTable().from_file("tmp/recur_serial");
  RecursiveTable recur_parallel = RecursiveTable().from_file("tmp/recur");
  EXPECT_TRUE(recur_serial.is_equal(recur_parallel));
}

TEST(BuildRecursiveTable, parallel_en_VERY_LONG) {
  const int num_proc = 3;

  argtype config_args = {{"side_length", "30,30,30"}, {"periodic", "false,false,false"}, {"particle_type", "t1:../plugin/aniso/particle/aniso_l10_patch.txt,t2:../particle/hard_sphere.txt"}, {"add_num_t1_particles", "1"}, {"add_num_t2_particles", "1"}, {"cutoff", str(std::pow(2, 1./6.))}, {"cutoffP_H", "3"}};
  int dimension = 2;
  std::string min_criteria_energy = "1.25", min_criteria = "0.25", ref_sig = "sigmaA_H", num_z = "5";
  if (false) {
  //if (true) {
    // 3d l - l interaction instead of default l-sphere
    config_args = {{"side_length", "50,50,50"}, {"periodic", "false,false,false"}, {"particle_type", "t1:../plugin/aniso/particle/aniso_l10_patch.txt"}, {"add_num_t1_particles", "2"}, {"cutoff", str(std::pow(2, 1./6.))}, {"cutoffP_P", "3"}};
    dimension = 5;
    min_criteria_energy = "2.5";
    min_criteria = "5";
    ref_sig = "sigmaA_A";
    num_z = "5";
  }

  // mayer sampling training data
  argtype mayer_config_args = config_args;
  mayer_config_args.insert({"group", "fixed,mobile"});
  mayer_config_args.insert({"fixed_particle_index", "0"});
  mayer_config_args.insert({"mobile_particle_index", "1"});
  //mayer_config_args.insert({"set_cutoff_min_to_sigma", "true"}); // for hard spheres

if (false) {
//if (true) {
  WARN("comment above to generate necessary files");
} else {

  MakeMonteCarlo({{
    {"Configuration", mayer_config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"RefPotential", {{"ref", "hs"}, {"Model", "HardSphere"}, {"sigma", "0"}, {ref_sig, "1"}}},
    {"ThermoParams", {{"beta", "1"}}},
    {"MayerSampling", {{"trials_per_cycle", "1e4"}, {"cycles_to_complete", "1e2"}}},
    {"TrialTranslate", {{"new_only", "true"}, {"ref", "hs"}, {"tunable_param", "1"}, {"group", "mobile"}}},
    {"TrialRotate", {{"new_only", "true"}, {"ref", "hs"}, {"tunable_param", "1"}, {"group", "mobile"}}},
    {"Tune", {{}}},
    {"Run", {{"until", "complete"}}},
    {"Remove", {{"name", "Tune"}}}, //,CriteriaWriter,Log,Movie
    {"CriteriaWriter", {{"trials_per_write", "1e4"}, {"output_file", "tmp/2recurb2.csv"}}},
    {"MayerSampling", {{"trials_per_cycle", "1e4"}, {"cycles_to_complete", "1e2"}, {"training_file", "tmp/2recur_training.csv"}}},
    {"Run", {{"until", "complete"}}},
  }}, true);

  // build base cutoff and contact
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage","base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2base_serial"}}},
  }}, true);
  RecursiveTable base_serial = RecursiveTable().from_file("tmp/2base_serial");
  for (int proc = 0; proc < num_proc; ++proc) {
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2base"+str(proc)}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
  }
  { // temporary testing
    RecursiveTable base_parallel0 = RecursiveTable().from_file("tmp/2base0");
    ASSERT(base_parallel0.cutoff2d_.size() > 0, "error");
    INFO(base_parallel0.cutoff2d_.size());
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/2base"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  RecursiveTable base_parallel = RecursiveTable().from_file("tmp/2base");
  EXPECT_TRUE(base_serial.is_equal(base_parallel));
  if (dimension == 2) {
    ASSERT(base_serial.cutoff2d_.size() > 0, "error");
    ASSERT(base_parallel.cutoff2d_.size() > 0, "error");
    EXPECT_NE(base_serial.contact2d_[0][0].data()[0][0],
              base_serial. cutoff2d_[0][0].data()[0][0]);
    EXPECT_NE(base_parallel.contact2d_[0][0].data()[0][0],
              base_parallel. cutoff2d_[0][0].data()[0][0]);
  } else if (dimension == 5) {
    ASSERT(base_serial.cutoff_.size() > 0, "error");
    ASSERT(base_parallel.cutoff_.size() > 0, "error");
    EXPECT_NE(base_serial.contact_[0][0].data()[0][0][0][0][0],
              base_serial. cutoff_[0][0].data()[0][0][0][0][0]);
    EXPECT_NE(base_parallel.contact_[0][0].data()[0][0][0][0][0],
              base_parallel. cutoff_[0][0].data()[0][0][0][0][0]);
  }

  // train contact and cutoff
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage", "train"}, {"mayer_training_file", "tmp/2recur_training.csv"}, {"base_table_file", "tmp/2base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2train_serial"}, {"min_criteria", min_criteria}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "train"}, {"mayer_training_file", "tmp/2recur_training.csv"}, {"base_table_file", "tmp/2base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2train"+str(proc)}, {"min_criteria", min_criteria}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "trained"}, {"input_prefix", "tmp/2train"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "trained"}, {"input_prefix", "tmp/2train"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}, {"input_suffix", "_cut"}}}}}, true);
  EXPECT_TRUE(compare_files("tmp/2train", "tmp/2train_serial"));
  EXPECT_TRUE(compare_files("tmp/2train_cut", "tmp/2train_serial_cut"));
  EXPECT_FALSE(compare_files("tmp/2train_cut", "tmp/2train"));

  // build nested contact table
  WARN("must somehow separate the trained file for contact with the trained file for cutoff");
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage", "build_contact"}, {"base_table_file", "tmp/2base"}, {"trained_file", "tmp/2train"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2recur_serial"}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    INFO("proc " << proc);
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "build_contact"}, {"base_table_file", "tmp/2base"}, {"trained_file", "tmp/2train"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2recur"+str(proc)}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
    INFO("combining");
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/2recur"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  //EXPECT_TRUE(compare_files("tmp/2recur", "tmp/2recur_serial"));

  RecursiveTable recur_serial = RecursiveTable().from_file("tmp/2recur_serial");
  RecursiveTable recur_parallel = RecursiveTable().from_file("tmp/2recur");
  EXPECT_TRUE(recur_serial.is_equal(recur_parallel));

  if (dimension == 5) {
    EXPECT_TRUE(recur_parallel.contact_.size() > 0);
  } else if (dimension == 2) {
    EXPECT_TRUE(recur_parallel.contact2d_.size() > 0);
  }
  // build nested cutoff table
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage", "build_cutoff"}, {"base_table_file", "tmp/2recur_serial"}, {"trained_file", "tmp/2train_serial_cut"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2recur_serial_cut"}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    INFO("proc " << proc);
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "build_cutoff"}, {"base_table_file", "tmp/2recur"}, {"trained_file", "tmp/2train_cut"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2recur_cut"+str(proc)}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
    INFO("combining");
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/2recur_cut"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);

  RecursiveTable recur_serial_cut = RecursiveTable().from_file("tmp/2recur_serial_cut");
  RecursiveTable recur_parallel_cut = RecursiveTable().from_file("tmp/2recur_cut");
  EXPECT_TRUE(recur_serial_cut.is_equal(recur_parallel_cut));

  if (dimension == 5) {
    EXPECT_TRUE(recur_parallel_cut.contact_.size() > 0);
    EXPECT_TRUE(recur_parallel_cut.cutoff_ .size() > 0);
  } else if (dimension == 2) {
    EXPECT_TRUE(recur_parallel_cut.contact2d_.size() > 0);
    EXPECT_TRUE(recur_parallel_cut.cutoff2d_ .size() > 0);
  }
}
  // now build energy base, using the recursive contact and cutoff
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage","base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2en_base_serial"}, {"energy_only", "true"}, {"num_z", num_z}, {"input_recursive_table", "tmp/2recur_serial_cut"}}},
  }}, true);
  RecursiveTable en_base_serial = RecursiveTable().from_file("tmp/2en_base_serial");
  for (int proc = 0; proc < num_proc; ++proc) {
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2en_base"+str(proc)}, {"energy_only", "true"}, {"num_z", num_z}, {"input_recursive_table", "tmp/2recur_cut"}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/2en_base"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  RecursiveTable en_base_parallel = RecursiveTable().from_file("tmp/2en_base");
  if (dimension == 5) {
    INFO("serial " << feasst_str(en_base_serial.energy_[0][0].data()[0][0][0][0][0]));
    INFO("parallel " << feasst_str(en_base_parallel.energy_[0][0].data()[0][0][0][0][0]));
  } else if (dimension == 2) {
    INFO("serial " << feasst_str(en_base_serial.energy3d_[0][0].data()[0][0]));
    INFO("parallel " << feasst_str(en_base_parallel.energy3d_[0][0].data()[0][0]));
    { // temporary testing
      for (int proc = 0; proc < num_proc; ++proc) {
        RecursiveTable ebp = RecursiveTable().from_file("tmp/2en_base"+str(proc));
        INFO(proc << " " << feasst_str(ebp.energy3d_[0][0].data()[0][0]));
      }
    }
  }
  ASSERT(en_base_serial.is_equal(en_base_parallel), "err");

  // train energy
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage", "train"}, {"mayer_training_file", "tmp/2recur_training.csv"}, {"base_table_file", "tmp/2en_base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2en_train_serial"}, {"num_z", num_z}, {"min_criteria", min_criteria}, {"min_criteria_energy", min_criteria_energy}, {"energy_only", "true"}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "train"}, {"mayer_training_file", "tmp/2recur_training.csv"}, {"base_table_file", "tmp/2en_base"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2en_train"+str(proc)}, {"num_z", num_z}, {"min_criteria", min_criteria}, {"min_criteria_energy", min_criteria_energy}, {"energy_only", "true"}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "trained"}, {"input_prefix", "tmp/2en_train"}, {"input_suffix", "_en"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  EXPECT_TRUE(compare_files("tmp/2en_train_en", "tmp/2en_train_serial_en"));

  // build recursive energy
  MakeMonteCarlo({{
    {"Configuration", config_args},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    {"BuildRecursiveTable", {{"stage", "build_energy"}, {"base_table_file", "tmp/2en_base_serial"}, {"trained_file", "tmp/2en_train_serial_en"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2en_recur_serial"}, {"num_z", num_z}, {"energy_only", "true"}}},
  }}, true);
  for (int proc = 0; proc < num_proc; ++proc) {
    INFO("proc " << proc);
    MakeMonteCarlo({{
      {"Configuration", config_args},
      {"Potential", {{"Model", "LennardJonesCutShift"}}},
      {"BuildRecursiveTable", {{"stage", "build_energy"}, {"base_table_file", "tmp/2en_base"}, {"trained_file", "tmp/2en_train_en"}, {"num_orientations_per_pi", "2"}, {"output_file", "tmp/2en_recur"+str(proc)}, {"num_z", num_z}, {"energy_only", "true"}, {"processor", str(proc)}, {"num_processors", str(num_proc)}}},
    }}, true);
    INFO("combining");
  }
  MakeMonteCarlo({{{"CombineParallelFiles", {{"type", "recursive_table"}, {"input_prefix", "tmp/2en_recur"}, {"num_processors", str(num_proc)}, {"cleanup", "false"}}}}}, true);
  //EXPECT_TRUE(compare_files("tmp/2recur", "tmp/2recur_serial"));
  RecursiveTable en_recur_serial = RecursiveTable().from_file("tmp/2en_recur_serial");
  RecursiveTable en_recur_parallel = RecursiveTable().from_file("tmp/2en_recur");
  EXPECT_TRUE(en_recur_serial.is_equal(en_recur_parallel));

  if (dimension == 2) {
    EXPECT_TRUE(en_recur_parallel.contact2d_[0][0].percent_nested() > 0);
    EXPECT_TRUE(en_recur_parallel.cutoff2d_ [0][0].percent_nested() > 0);
    EXPECT_TRUE(en_recur_parallel.energy3d_ [0][0].percent_nested() > 0);
  } else if (dimension == 5) {
    EXPECT_TRUE(en_recur_parallel.contact_[0][0].percent_nested() > 0);
    EXPECT_TRUE(en_recur_parallel.cutoff_ [0][0].percent_nested() > 0);
    EXPECT_TRUE(en_recur_parallel.energy_ [0][0].percent_nested() > 0);
  }
}

}  // namespace feasst
