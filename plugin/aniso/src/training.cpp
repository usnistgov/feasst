#include <iostream>
#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "utils/include/io.h"
#include "utils/include/utils.h" // resize
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/position.h"
#include "math/include/matrix.h"
#include "math/include/euler.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "configuration/include/model_params.h"
#include "system/include/system.h"
#include "system/include/hard_sphere.h"
#include "system/include/visit_model.h"
#include "system/include/visit_model_inner.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/training.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

Training::Training(argtype * args) {
  class_name_ = "Training";
  training_file_ = str("training_file", args);
  output_file_ = str("output_file", args, "");
}
Training::Training(argtype args) : Training(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(Training, argtype({{"training_file", "tmp"}}));

Training::Training(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6274, "mismatch version: " << version);
}

void Training::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(6274, ostr);
  WARN("file names are not serialized");
}

void Training::run(MonteCarlo * mc) {
  const System& sys = mc->system();
  ASSERT(sys.num_configurations() == 1, "Assumes 1 Configuration");
  const Configuration& config = sys.configuration();
  ASSERT(config.dimension() == 3, "Assumes 3D");
  ASSERT(config.num_particles() == 2, "Assumes 2 particles");
  std::ifstream file(training_file_.c_str());
  std::ofstream output;
  if (!output_file_.empty()) {
    output.open(output_file_.c_str());
  }
  ASSERT(file.good(), "cannot find file " << training_file_.c_str());
  std::string line = "init";
  bool last_line = false;
  int num = 0;
  int match = 0;
  int hard = 0;
  int soft = 0;
  Position cart; //, rel; //, pbc, pos1, sph;
  cart.set_to_origin(3);
  WARN("types hardcoded to be the same.");
  int type1 = 0, type2 = 0;
  //rel.set_to_origin(3);
  //pbc.set_to_origin(3);
  //pos1.set_to_origin(3);
  //sph.set_to_origin(3);
  double s1, s2, e1, e2, e3;
  Euler euler;
  std::stringstream ss;
  mc->system().potential(0).visit_model().inner().serialize(ss);
  auto inner = std::make_unique<VisitModelInnerTable>(ss);
  const Table5D& contact = *config.table5d()[0][0];
  struct stat {
    int num;
    int match;
    int soft;
    int hard;
  };
  std::vector<std::vector<std::vector<std::vector<std::vector<stat> > > > > dat;
  resize(contact.num(0), contact.num(1), contact.num(2), contact.num(3), contact.num(4), &dat);
  //const std::vector<std::vector<std::shared_ptr<Table5D> > >& contact = config.table5d();
  INFO("contact size " << contact.num0() << " x " << contact.num1() << " x " << contact.num2() << " x " << contact.num3() << " x " << contact.num4());
  //std::vector<std::vector<std::shared_ptr<Table5D> > > * contact = mc->get_system()->get_configuration()->get_table5d();
  //HardSphere model;
  //RotationMatrix rot;
  //const double cutoff = mc->system().configuration().model_params().select("cutoff").value(0);
  WARN("why flip particles in config update positions?");
  while (!line.empty() && !last_line) {
    if (file.eof()) last_line = true;
    std::getline(file, line);
    if (!line.empty()) {
      ++num;
      std::vector<std::string> data = split(line, ',');
      const double srho = str_to_double(data[0]);
      const double stheta = str_to_double(data[1]);
      const double sphi = str_to_double(data[2]);
      cart.set_from_spherical(srho, stheta, sphi);
      DEBUG("cart: " << cart.str());
      const double ephi = str_to_double(data[3]);
      const double etheta = str_to_double(data[4]);
      const double epsi = str_to_double(data[5]);
      const double exp = str_to_double(data[6]);
      //INFO("computing");
      mc->get_system()->get_configuration()->update_positions(
        {{cart.coord(0), cart.coord(1), cart.coord(2)}, {0, 0, 0}},
        {{ephi, etheta, epsi}, {0, 0, 0}});
//        {{0, 0, 0}, {cart.coord(0), cart.coord(1), cart.coord(2)}},
//        {{0, 0, 0}, {ephi, etheta, epsi}});
      const double en = mc->get_system()->energy();
      //INFO("recomputing");
      inner->compute_scaled_coords(0, 1,
        config.particle(0).site(0), config.particle(1).site(0), const_cast<const Configuration *>(&config),
        &type1, &type2, &cart, &s1, &s2, &e1, &e2, &e3);
      //INFO("s1 " << s1 << " s2 " << s2 << " e1 " << e1 << " e2 " << e2 << " e3 " << e3);
      //INFO("tab " << contact.value_to_lowest_bin(0, s1));
      //INFO("tab " << contact.value_to_nearest_bin(0, s1));
      const int in0 = contact.value_to_lowest_bin(0, s1);
      const int in1 = contact.value_to_lowest_bin(1, s2);
      const int in2 = contact.value_to_lowest_bin(2, e1);
      const int in3 = contact.value_to_lowest_bin(3, e2);
      const int in4 = contact.value_to_lowest_bin(4, e3);
      const double en2 = inner->compute_aniso(type1, type2, srho*srho, s1, s2, e1, e2, e3, config);
      stat * st = &dat[in0][in1][in2][in3][in4];
      ASSERT(std::abs(en - en2) < NEAR_ZERO, "en: " << en << " en2: " << en2);
      st->num += 1;
      if ( (en == 0 && exp == 0) ||
           (en > 0 && exp > 0) ) {
        ++match;
        st->match += 1;
      } else if (en > exp) {
        ++hard;
        st->hard += 1;
      } else if (en <= exp) {
        ++soft;
        st->soft += 1;
      } else {
        FATAL("unreachable");
      }
      //DEBUG(config.particle(1).site(0).position().str() << " " << feasst_str(data) << " en:" << en << " exp:" << exp << " match:" << match);
    }
  }
  INFO("percent correct:" << static_cast<double>(match)/num
    << " hard:" << static_cast<double>(hard)/num
    << " soft:" << static_cast<double>(soft)/num
  );
  if (output.good()) {
    for (const auto& st1 : dat) {
    for (const auto& st2 : st1) {
    for (const auto& st3 : st2) {
    for (const auto& st4 : st3) {
    for (const auto& st5 : st4) {
      if (st5.num != 0) {
        output << st5.num << " " << static_cast<double>(st5.match)/st5.num
               << " " << static_cast<double>(st5.hard)/st5.num
               << " " << static_cast<double>(st5.soft)/st5.num
               << std::endl;
      }
    }}}}}
  }
}

}  // namespace feasst
