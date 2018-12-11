#include <iostream>
#include "feasst.h"

int main() {
  feasst::Position pos;
  pos.set_to_origin_3D();
  std::vector<double> x = {3.5, 796.4, -45.4};
  pos.set_vector(x);
  std::vector<double> x2 = pos.coord();
  std::cout << pos.coord(0) << std::endl;
  feasst::Configuration config;
}
