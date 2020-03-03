
#include <fstream>
#include "xtc/include/file_xtc.h"

namespace feasst {

int FileXTC::load(XDRFILE * file,
                  Configuration * config,
                  std::string file_name) {
  const int dimension = 3;
  int endXTC = 0;
//  XDRFILE * xdr_file = xdrfile_open(file_name.c_str(), "r");
  // check atoms
  int natoms_xtc;
  int result_xtc;
  char * fn_xtc = new char[1];
  if (file_name.empty()) {
    natoms_xtc = config->num_sites();
  } else {
    delete [] fn_xtc;
    fn_xtc = new char[file_name.size() + 1];
    std::copy(file_name.begin(), file_name.end(), fn_xtc);
    fn_xtc[file_name.size()] = '\0';
    result_xtc = read_xtc_natoms(fn_xtc, &natoms_xtc);
    ASSERT(exdrOK == result_xtc,
      "cannot read natoms, " << natoms_xtc << " from xtc file: " << file_name);
    ASSERT(config->num_sites() == natoms_xtc, "number of atoms(" << config->num_sites()
      << ") does not match number of atoms in xtc, " << natoms_xtc
      << " from xtc file: " << file_name);
  }

  // set coordinates and box length
  int step_xtc;
  float time_xtc;
  matrix box_xtc;
  rvec *x_xtc;
  x_xtc = reinterpret_cast<rvec *>(calloc(natoms_xtc, sizeof(x_xtc[0])));
  float prec_xtc = 1000.0;
  result_xtc = read_xtc(file, natoms_xtc, &step_xtc, &time_xtc,
                        box_xtc, x_xtc, &prec_xtc);
  if (result_xtc != 0) {
    INFO("reached the end of XTC file " << file_name);
    endXTC = 1;
  }
  // set the domain
  DEBUG("box: " << box_xtc[0][0] << " "
                << box_xtc[1][1] << " "
                << box_xtc[2][2] << " "
                << box_xtc[0][1]);
  config->set_side_lengths(
    Position().set_vector({box_xtc[0][0], box_xtc[1][1], box_xtc[2][2]}));

  // set the coordinates
  std::vector<std::vector<double> > coords;
  coords.resize(natoms_xtc, std::vector<double>(dimension));
  for (int i = 0; i < natoms_xtc; ++i) {
    for (int dim = 0; dim < dimension; ++dim) {
      coords[i][dim] = x_xtc[i][dim];
    }
  }
  config->update_positions(coords);
  free(x_xtc);
  if (!file_name.empty()) {
    delete [] fn_xtc;
  }
  return endXTC;
}

void FileXTC::write(XDRFILE * file, const Configuration& config) {
  int natoms_xtc = config.num_sites();
  matrix box_xtc;
  box_xtc[0][0] = config.domain()->side_length(0);
  box_xtc[0][1] = 0;
  box_xtc[0][2] = 0;
  box_xtc[1][0] = 0;
  box_xtc[1][1] = config.domain()->side_length(1);
  box_xtc[1][2] = 0;
  box_xtc[2][0] = 0;
  box_xtc[2][1] = 0;
  box_xtc[2][2] = config.domain()->side_length(2);
  rvec *x_xtc;
  x_xtc = reinterpret_cast<rvec *>(calloc(natoms_xtc, sizeof(x_xtc[0])));
  float prec_xtc = 1000.0;
  int atom = 0;
  for (const Particle& part : config.particles().particles()) {
    for (const Site& site : part.sites()) {
  //for (int i = 0; i < natoms_xtc; ++i) {
      for (int dim = 0; dim < config.dimension(); ++dim) {
        x_xtc[atom][dim] = site.position().coord(dim);
      }
      ++atom;
    }
  }
  if (write_xtc(file, natoms_xtc, 0, 0, box_xtc, x_xtc, prec_xtc) != 0) {
    ERROR("error writing xtc file");
  }
  free(x_xtc);
}

}  // namespace feasst
