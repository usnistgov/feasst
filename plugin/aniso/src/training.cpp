#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "utils/include/io.h"
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
}

void Training::run(MonteCarlo * mc) {
  const System& sys = mc->system();
  ASSERT(sys.num_configurations() == 1, "Assumes 1 Configuration");
  const Configuration& config = sys.configuration();
  ASSERT(config.dimension() == 3, "Assumes 3D");
  ASSERT(config.num_particles() == 2, "Assumes 2 particles");
  std::ifstream file(training_file_.c_str());
  ASSERT(file.good(), "cannot find file " << training_file_.c_str());
  std::string line = "init";
  bool last_line = false;
  int num = 0;
  int match = 0;
  int hard = 0;
  int soft = 0;
  Position cart, rel, pbc, pos1, sph;
  cart.set_to_origin(3);
  rel.set_to_origin(3);
  pbc.set_to_origin(3);
  pos1.set_to_origin(3);
  sph.set_to_origin(3);
  Euler euler;
  //VisitModelInner * tab = mc->get_system()->get_potential(0)->get_visit_model_()->get_inner_();
  std::stringstream ss;
  mc->system().potential(0).visit_model().inner().serialize(ss);
  auto tab = std::make_unique<VisitModelInnerTable>(ss);
  HardSphere model;
  RotationMatrix rot;
  const double cutoff = mc->system().configuration().model_params().select("cutoff").value(0);
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
      mc->get_system()->get_configuration()->update_positions(
        {{0, 0, 0}, {cart.coord(0), cart.coord(1), cart.coord(2)}},
        {{0, 0, 0}, {ephi, etheta, epsi}});
      //const double en = mc->initialize_system(0);
      const double en = mc->get_system()->energy();
      
//      // flip and scale the coordinates
//      euler.set(ephi, etheta, epsi);
//      euler.compute_rotation_matrix(&rot);
//      rot.transpose();
//      rot.multiply(cart, &pos1);
//      pos1.spherical(&sph);
//      euler.set(rot);
//      double sctheta = sph.coord(1)/PI;
//      //double sctheta = stheta/PI;
//      if (sctheta < 0) sctheta += 1;
//      const double scphi = sph.coord(2)/PI;
//      //const double scphi = sphi/PI;
//      const double ecphi = euler.phi()/2/PI + 0.5;
//      //const double ecphi = ephi/2/PI + 0.5;
//      const double ectheta = euler.theta()/PI;
//      //const double ectheta = etheta/PI;
//      const double ecpsi = euler.psi()/2/PI + 0.5;
//      //const double ecpsi = epsi/2/PI + 0.5;
//      DEBUG(srho*srho << " " << cutoff*cutoff);
//      const double en2 = tab->compute_aniso(0, 0, srho*srho, sctheta, scphi, ecphi, ectheta, ecpsi, config);
//      //const double en2 = tab->compute_aniso(0, 0, srho*srho, stheta, sphi, ephi, etheta, epsi, config);
//      //tab->compute(0, 0, 1, 0, mc->get_system()->get_configuration(), config.model_params(), &model, 0, &rel, &pbc, 1);
//      //const double en2 = tab->energy();
//      ASSERT(std::abs(en - en2) < NEAR_ZERO, "en: " << en << " en2: " << en2);
      if ( (en == 0 && exp == 0) ||
           (en > 0 && exp > 0) ) {
        ++match;
      } else if (en > exp) {
        ++hard;
      } else if (en <= exp) {
        ++soft;
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
}

}  // namespace feasst
