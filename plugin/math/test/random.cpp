#include <cmath>
#include <vector>
#include <fstream>
#include "utils/test/utils.h"
#include "math/include/random_modulo.h"
#include "math/include/random_mt19937.h"
#include "math/include/histogram.h"
#include "math/include/accumulator.h"
#include "math/include/matrix.h"

namespace feasst {

std::vector<std::shared_ptr<Random> > gens = {MakeRandomMT19937(), MakeRandomModulo()};
//std::vector<std::shared_ptr<Random> > gens = {MakeRandomModulo()};

TEST(Random, uniform) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    double sum = 0;
    const int num = 1e5;
    for (int i = 0; i < num; ++i) {
      const double uni = random->uniform();
      //INFO("uni " << uni);
      sum += uni;
    }
    EXPECT_NEAR(sum, double(num)/2., sqrt(double(num)));
  }
}

TEST(Random, uniform_int) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    double sum = 0;
    const int num = 1e5;
    for (int i = 0; i < num; ++i) {
      sum += random->uniform(0, 10);
    }
    EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
  }
}

TEST(Random, alpha_numeric) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    const int size = 10;
    std::string unique = random->alpha_numeric(size);
    // INFO("unique alpha numeric: " << unique);
    EXPECT_EQ(unique.size(), size);
    EXPECT_NE(unique, random->alpha_numeric(size));
  }
}

TEST(Random, element) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    std::vector<int> vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double sum = 0;
    const int num = 1e5;
    for (int i = 0; i < num; ++i) {
      sum += random->const_element(vec);
    }
    EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
  }
}

TEST(Random, unit_sphere) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    Position position;
    position.set_vector({0., 0.});
    random->unit_sphere_surface(&position);
    EXPECT_NEAR(position.distance(), 1., NEAR_ZERO);
    position.set_vector({0., 0., 0.});
    random->unit_sphere_surface(&position);
    EXPECT_NEAR(position.distance(), 1., NEAR_ZERO);

  //  // visualize
  //  for (int point = 0; point < 1e3; ++point) {
  //    random.unit_sphere_surface(&position);
  //    std::cout << position.str() << std::endl;
  //  }
  }
}

TEST(Random, spherical_shell) {
  const double upper = 3, lower = 1.;
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    Position position;

    // 3D
    position.set_vector({0., 0., 0.});
    for (int point = 0; point < 1e3; ++point) {
      random->position_in_spherical_shell(lower, upper, &position);
      EXPECT_LE(lower, position.distance());
      EXPECT_GE(upper, position.distance());
      // visualize
      // std::cout << position.str() << std::endl;
    }

    // 2D
    position.set_vector({0., 0.});
    for (int point = 0; point < 1e3; ++point) {
      random->position_in_spherical_shell(lower, upper, &position);
      EXPECT_LE(lower, position.distance());
      EXPECT_GE(upper, position.distance());
      // visualize
      // std::cout << position.str() << std::endl;
    }
  }
}

TEST(Random, index_from_cumulative_probability) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    std::vector<double> cpdf;
    const int ncpdf = 10, num = 100;
    for (int i = 0; i < ncpdf; ++i) cpdf.push_back((i+1)/double(ncpdf));
    std::vector<double> cpdfran(ncpdf);
    for (int i = 0; i < num; ++i) {
      const int j = random->index_from_cumulative_probability(cpdf);
      ++cpdfran[j];
    }
    for (int i = 0; i < ncpdf; ++i) EXPECT_NEAR(cpdfran[i]/double(num), ncpdf/double(num), 0.2);
  }
}

TEST(RandomModulo, compiler_independent) {
  auto random = MakeRandomModulo({{"seed", "1346867550"}});
//  EXPECT_NEAR(random.uniform(), 0, NEAR_ZERO);
//  EXPECT_NEAR(random.uniform(), 0, NEAR_ZERO);
//  EXPECT_NEAR(random.uniform(), 0, NEAR_ZERO);
  EXPECT_EQ(random->uniform(), 0.082789837886947132);
  EXPECT_EQ(random->uniform(), 0.44880536592044185);
  EXPECT_EQ(random->uniform(), 0.071785024866361652);
}

TEST(RandomModulo, serialize) {
  RandomModulo random;
  random.seed_by_time();
  random.set_cache_to_load(true);
  RandomModulo random2 = test_serialize(random);
  const double next = random.uniform();
  EXPECT_EQ(next, random2.uniform());
  RandomModulo random3;
  random3.set_cache_to_unload(random2);
  EXPECT_EQ(next, random3.uniform());
  random3.set_cache_to_load(false);
  EXPECT_NE(random.uniform(), random3.uniform());
  TRY(
    random3.set_cache_to_unload(random2);
    random3.uniform();
    random3.uniform();
    CATCH_PHRASE("can not unload if nothing stored");
  );
}

TEST(RandomMT19937, serialize) {
  RandomMT19937 random;
  random.seed_by_time();
  random.set_cache_to_load(true);
  RandomMT19937 random2 = test_serialize(random);
  const double next = random.uniform();
  EXPECT_EQ(next, random2.uniform());
  RandomMT19937 random3;
  random3.set_cache_to_unload(random2);
  EXPECT_EQ(next, random3.uniform());
  random3.set_cache_to_load(false);
  EXPECT_NE(random.uniform(), random3.uniform());
  TRY(
    random3.set_cache_to_unload(random2);
    random3.uniform();
    random3.uniform();
    CATCH_PHRASE("can not unload if nothing stored");
  );
}

TEST(Random, standard_normal) {
  for (std::shared_ptr<Random> random : gens) {
    random->seed_by_time();
    double sum = 0., sumsq = 0.;
    const int num = 1e5;
    const double num_dble = static_cast<double>(num);
    for (int i = 0; i < num; ++i) {
      const double ran = random->standard_normal();
      sum += ran;
      sumsq += ran*ran;
    }
    const double av = sum/num_dble;
    const double stdev = std::sqrt((sumsq/num_dble-av*av));
    EXPECT_NEAR(sum, 0., 4*std::sqrt(num_dble));
    EXPECT_NEAR(stdev, 1., 0.01);
  }
}

TEST(Random, rotation_tunable) {
  auto ran = MakeRandomMT19937();
  TRY(
    RotationMatrix rotmat = ran->rotation(3, -0.5);
    CATCH_PHRASE("must be -1 or > 0");
  );
  RotationMatrix rotmat = ran->rotation(3, 0);
}

TEST(Random, rotation) {
  std::ofstream file("tmp/tmp.txt");
  auto ran = MakeRandomMT19937();
  const Position origin({0., 0., 0.});
  //const double tunable = 1;
  const double tunable = -1;
  for (int trial = 0; trial < 1e4; ++trial) {
    RotationMatrix rotmat = ran->rotation(3, tunable);
    Position unit({1., 0., 0.});
    rotmat.rotate(origin, &unit);
    file << unit.coord(0) << " " << unit.coord(1) << " " << unit.coord(2) << std::endl;
  }
}

//TEST(Random, bond_length) {
//  for (std::shared_ptr<Random> random : gens) {
//    random->seed_by_time();
//    double spring_constant = 100.;
//    const double equilibrum_length = 3.;
//    double sum1 = 0., sumsq1 = 0., sum2 = 0., sumsq2 = 0.;
//    const int num = 5e3;
//    const double num_dble = static_cast<double>(num);
//  //  Histogram hist1, hist2;
//  //  hist1.set_width_center(0.01, equilibrum_length);
//  //  hist2.set_width_center(0.01, equilibrum_length);
//    //({{"width", "0.1"}, {"max", "6"}});;
//    for (int i = 0; i < num; ++i) {
//      const double ran1 = random->bond_length(equilibrum_length, 2*equilibrum_length, spring_constant, 2, 3);
//      const double ran2 = random->harmonic_bond_length(equilibrum_length, spring_constant, 3);
//  //    hist1.add(ran1);
//  //    hist2.add(ran2);
//      sum1 += ran1; sumsq1 += ran1*ran1; sum2 += ran2; sumsq2 += ran2*ran2;
//    }
//    const double av1 = sum1/num_dble;
//    const double av2 = sum2/num_dble;
//    const double stdev1 = std::sqrt((sumsq1/num_dble-av1*av1));
//    const double stdev2 = std::sqrt((sumsq2/num_dble-av2*av2));
//    const double stdev_exp = 1./std::sqrt(2*spring_constant);
//    EXPECT_NEAR(av1, equilibrum_length, 8e-2);
//    EXPECT_NEAR(av2, equilibrum_length, 8e-2);
//    EXPECT_NEAR(stdev1, stdev_exp, 5e-2);
//    EXPECT_NEAR(stdev2, stdev_exp, 5e-2);
//  //  for (int bin = 0; bin < hist.size(); ++bin) {
//  //    std::cout << hist.center_of_bin(bin) << " " << hist.histogram()[bin] << std::endl;
//  //  }
//  }
//}

//TEST(Random, bond_angle) {
//  for (std::shared_ptr<Random> random : gens) {
//    random->seed_by_time();
//    double sum1 = 0., sumsq1 = 0., sum2 = 0., sumsq2 = 0.;
//    const int num = 2e3;
//    const double num_dble = static_cast<double>(num);
//    Position point(3);
//    for (int i = 0; i < num; ++i) {
//      const double ran1 = random->bond_angle(0, 0, 2, 3);
//      random->position_in_spherical_shell(0, 1., &point);
//      const double ran2 = point.spherical().coord(2);
//      sum1 += ran1; sumsq1 += ran1*ran1; sum2 += ran2; sumsq2 += ran2*ran2;
//    }
//    const double av1 = sum1/num_dble;
//    const double av2 = sum2/num_dble;
//    const double stdev1 = std::sqrt((sumsq1/num_dble-av1*av1));
//    const double stdev2 = std::sqrt((sumsq2/num_dble-av2*av2));
//    const double stdev_exp = 0.675;
//    EXPECT_NEAR(av1, PI/2, 8e-2);
//    EXPECT_NEAR(av2, PI/2, 8e-2);
//    EXPECT_NEAR(stdev1, stdev_exp, 5e-2);
//    EXPECT_NEAR(stdev2, stdev_exp, 5e-2);
//  }
//}

// Compute the average bond length for a potential U(L)=600(L-1)^2
TEST(Random, harmonic_bond_brute) {
  const double k = 600;
  const double l0 = 1;
  auto ran = MakeRandomMT19937();
  double old_dist = l0;
  Position xn({{"x", feasst::str(old_dist)}, {"y", "0"}, {"z", "0"}});
  Position old_x = xn;
  Accumulator length;
  double old_en = 0.;
  const double max_move = 0.1;
  int accepted = 0;
  const int num_attempts = 3e5;
  //const int num_attempts = 1e7;
  for (int i = 0; i < num_attempts; ++i) {
    //ran->position_in_spherical_shell(0.5, 1.5, &xn);
    for (int dim = 0; dim < xn.dimension(); ++dim) {
      xn.set_coord(dim, old_x.coord(dim) + max_move*(2*ran->uniform() - 1));
    }
    const double dist = xn.distance(), dx = dist - l0;
    const double en = k*dx*dx;
    const double delta_en = en - old_en;
    if (ran->uniform() < std::exp(-delta_en)) {
      old_dist = dist;
      old_en = en;
      old_x = xn;
      ++accepted;
    }
    length.accumulate(old_dist);
  }
  DEBUG("acceptance " << static_cast<double>(accepted)/num_attempts);
  DEBUG(length.str());
//  for (int op = 0; op < length.max_block_operations(); ++op) {
//    if (length.block_averages()[op]) {
//      INFO(length.block_size()[op] << " numb " << length.block_averages()[op]->num_values() << " " << length.block_averages()[op]->average() << " " << length.block_stdev(op) << " " << length.block_std_of_std(op));
//    }
//  }
  EXPECT_NEAR(1.+1./k, length.average(), 15*length.block_stdev());
}

// Compute the average bond angle for a potential U(theta)=600(theta-PI/2)^2
TEST(Random, harmonic_angle_brute) {
  const double k = 600;
  const double theta0 = PI/2;
  double old_theta = theta0;
  auto ran = MakeRandomMT19937();
  Position xa1({{"x", "0"}, {"y", "1"}, {"z", "0"}});
  Position xa0({{"x", "0"}, {"y", "0"}, {"z", "0"}});
  Position xn({{"x", "1"}, {"y", "0"}, {"z", "0"}});
  Position old_x = xn;
  Accumulator angle;
  double old_en = 0.;
  const double max_move = 1;
  int accepted = 0;
  const int num_attempts = 3e5;
  //const int num_attempts = 1e7;
  for (int i = 0; i < num_attempts; ++i) {
    //ran->position_in_spherical_shell(0.5, 1.5, &xn);
    for (int dim = 0; dim < xn.dimension(); ++dim) {
      xn.set_coord(dim, old_x.coord(dim) + max_move*(2*ran->uniform() - 1));
    }
    const double theta = xa0.vertex_angle_radians(xn, xa1);
    const double dtheta = theta - theta0;
    const double en = k*dtheta*dtheta;
    const double delta_en = en - old_en;
    if (ran->uniform() < std::exp(-delta_en)) {
      old_theta = theta;
      old_en = en;
      old_x = xn;
      ++accepted;
    }
    angle.accumulate(old_theta);
  }
  DEBUG("acceptance " << static_cast<double>(accepted)/num_attempts);
  DEBUG(angle.str() << " PI/2 " << PI/2);
  EXPECT_NEAR(PI/2, angle.average(), 20*angle.block_stdev());
}

// Compute Henry coefficient of 3D HS dimer in slit pore.
// dimer bond length = sigma
// slit width = 6, box length = 10
// W-3/L is probably to insert inside slit where the second site is OK no matter what (1.5 sigma from either wall). 2/L is probably to insert first bead where second bead matters. Spherical cap surface area varies linearly with cap height, so its average acceptance probability of 1 or 1/2 which gives 3/4 factor
TEST(Random, henry_dimer_slit) {
  const double L = 10, W = 6, sigma = 1;
  Position site0;
  site0.set_to_origin(3);
  Position site1 = site0;
  Position quat;
  quat.set_to_origin(4);
  //bond = site0;
  const int dimen = 3;
  const int num = 3e5;
  auto ran = MakeRandomMT19937();
  Accumulator H;
  RotationMatrix R;
  for (int trial = 0; trial < num; ++trial) {
    bool inside = true;
    for (int dim = 0; dim < dimen; ++dim) {
      const double x= ran->uniform_real(-L/2, L/2);
      site0.set_coord(dim, x);
      if ((dim == 0) && (x > (W-sigma)/2 || x < -(W-sigma)/2)) inside = false;
    }
    if (inside) {
      //ran->unit_sphere_surface(&bond);
      //site1 = site0;
      //site1.add(bond);
      //const double ang = ran->uniform(-360, 360);
      //const double ang = ran->uniform(-180, 180);
      //ran->rotation(3, &bond, &R, ang);
      ran->rotation(3, &quat, &R);
      site1 = site0;
      site1.add_to_coord(0, 1);
      R.rotate(site0, &site1);

      const double x = site1.coord(0);
      if (x > (W-sigma)/2 || x < -(W-sigma)/2) inside = false;
    }
    if (inside) {
      H.accumulate(1);
      DEBUG("inside " << site0.str() << " " << site1.str());
    } else {
      H.accumulate(0);
      DEBUG("outside " << site0.str() << " " << site1.str());
    }
  }
  DEBUG(H.str());
  EXPECT_NEAR((W-3*sigma)/L+(2*sigma/L)*(3./4.), H.average(), 10*H.block_stdev());
}

// In this example, a random walk is performed in one dimension with periodic boundaries.
// The step size is chosen randomly uniform in [-1, 1], and the repeating cell is 25.
TEST(Accumulator, block_average_LONG) {
  const double half_pbc_length = 25/2.;
  const double step_size = 1;
  const int num_steps = 1e6;
  auto av_position = MakeAccumulator({{"num_moments", "5"}, {"max_block_operations", "20"}});
  auto av_position2 = MakeAccumulator({{"num_moments", "5"}, {"max_block_operations", "6"}});
  auto av_pos_sq = MakeAccumulator({{"num_moments", "5"}, {"max_block_operations", "20"}});
  DEBUG(feasst_str(av_position->block_size()));
  EXPECT_EQ(av_position->num_moments(), 5);
  //auto rng = MakeRandomMT19937({{"seed", "123"}});
  auto rng = MakeRandomMT19937();
  std::vector<double> pos;
  double position = rng->uniform_real(-half_pbc_length, half_pbc_length);
  for (int step = 0; step < num_steps; ++step) {
    position += step_size*rng->uniform_real(-1, 1);
    if (position < -half_pbc_length) {
      position += 2*half_pbc_length;
    } else if (position > half_pbc_length) {
      position -= 2*half_pbc_length;
    }
    av_position->accumulate(position);
    av_position2->accumulate(position);
    av_pos_sq->accumulate(std::sqrt(position*position));
    pos.push_back(position);
  }

  const double av = av_position->average();
  EXPECT_NEAR(av, std::accumulate(pos.begin(), pos.end(), 0.)/num_steps, 1e-8);
  //INFO(av_position->std());
  for (int mo = 2; mo <= 4; ++mo) {
    double mom = 0.;
    for (const double& p : pos) {
      mom += std::pow(p - av, mo);
    }
    mom /= num_steps;
    //INFO(mo << " " << mom << " " << av_position->central_moment(mo));
    EXPECT_NEAR(mom, av_position->central_moment(mo), 1e-8);
  }

  for (std::shared_ptr<Accumulator> acc : {av_position, av_position2}) {
    for (int op = 0; op < acc->max_block_operations(); ++op) {
      if (acc->block_averages()[op]) {
        INFO(acc->block_size()[op] << " numb " << acc->block_averages()[op]->num_values() << " " << acc->block_averages()[op]->average() << " " << acc->block_stdev(op) << " " << acc->block_std_of_std(op));
      }
    }
  }
//  // expect 1/3 L^2
//  for (int op = 0; op < av_position->max_block_operations(); ++op) {
//    INFO(op << " " << av_pos_sq->block_averages()[op]->average() << " " << av_pos_sq->block_stdev(op) << " " << av_pos_sq->block_std_of_std(op));
//  }
//  // consider an accumulator that is transformed to the final result. Blocking size can affect it
  INFO(av_position2->block_stdev());
  EXPECT_NEAR(av_position2->block_stdev(), 0.08, 0.05);
  EXPECT_NEAR(av_position->block_std_of_std(11), 0.05, 0.05);
}

}  // namespace feasst
