#include <iostream>
#include <omp.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <ctime>

#define DEBUG_SERIAL_MODE_5324634

// Simulation parameters
constexpr double cutoff = 3.;
constexpr int dimen = 3;
constexpr double max_move = 0.185;
//constexpr int num_trials = 2e4;
//constexpr int num_trials = 2e3; // same
//constexpr int num_trials = 20;
constexpr int num_trials = 1e6;
constexpr double beta = 1.2;
constexpr double box_length = 8;
constexpr int num_particles = 50;
constexpr bool verbose = false;
//constexpr bool verbose = true;

// Read an xyz file
std::vector<std::vector<double> > read_xyz(const std::string file_name = "../plugin/monte_carlo/test/data/bench.xyz") {
  std::vector<std::vector<double> > xyz(num_particles, std::vector<double>(dimen, 0));
  std::ifstream xyz_file(file_name);
  std::string line;
  std::getline(xyz_file, line);
  std::getline(xyz_file, line);
  for (std::vector<double>& pos : xyz) {
    std::getline(xyz_file, line);
    std::stringstream iss(line);
    std::string tmp;
    iss >> tmp;
    for (double& x : pos) {
      iss >> x;
    }
  }
  return xyz;
}

// Return the energy of a particle at a given position.
double en_part(const int ipart, const double xi, const double yi, const double zi,
  const std::vector<std::vector<double> >& xyz) {
  double en = 0.;
  double dx, dy, dz;
  for (int jpart = 0; jpart < num_particles; ++jpart) {
    if (ipart != jpart) {
      dx = xyz[jpart][0] - xi;
      if (dx > +0.5*box_length) dx -= box_length;
      if (dx < -0.5*box_length) dx += box_length;
      dy = xyz[jpart][1] - yi;
      if (dy > +0.5*box_length) dy -= box_length;
      if (dy < -0.5*box_length) dy += box_length;
      dz = xyz[jpart][2] - zi;
      if (dz > +0.5*box_length) dz -= box_length;
      if (dz < -0.5*box_length) dz += box_length;
      const double rsq = dx*dx + dy*dy + dz*dz;
      if (rsq < cutoff*cutoff) {
        const double rsq6inv = 1./(rsq*rsq*rsq);
        en += 4*rsq6inv*(rsq6inv - 1.);
      }
    }
  }
  return en;
}

// Return the energy of all particles.
double en_all(const std::vector<std::vector<double> >& xyz) {
  double en = 0.;
  for (int ipart = 0; ipart < num_particles; ++ipart) {
    const std::vector<double>& x = xyz[ipart];
    en += en_part(ipart, x[0], x[1], x[2], xyz);
  }
  return 0.5*en;
}

class RandomBenchmark {
 public:
  RandomBenchmark() { seed(); }
  void seed() {
    //srand(1234);  // seed for reproduction
    generator = std::mt19937(rand());
  }
  std::uniform_real_distribution<double> ran_uniform;
  std::mt19937 generator;
  double uniform() { return ran_uniform(generator); }
};

// Conduct a serial Metropolis Monte Carlo simulation.
void serial(RandomBenchmark * random) {
  std::vector<std::vector<double> > xyz = read_xyz();
  double energy = en_all(xyz);
  std::cout << "en all: " << energy << std::endl;
  double energy_sum = 0;
  int num_accepted = 0;
  for (int trial = 0; trial < num_trials; ++trial) {
    const int trial_particle = static_cast<int>(random->uniform()*num_particles);
    if (verbose) std::cout << "p " << trial_particle << std::endl;
    double xi = xyz[trial_particle][0];
    double yi = xyz[trial_particle][1];
    double zi = xyz[trial_particle][2];
    const double en_old = en_part(trial_particle, xi, yi, zi, xyz);
//    std::vector<double> xyz_old = xyz[trial_particle];
    double ran = random->uniform();
    xi += max_move*(ran - 0.5);
    if (verbose) std::cout << "ran " << ran <<std::endl;
    yi += max_move*(random->uniform() - 0.5);
    zi += max_move*(random->uniform() - 0.5);
    const double delta_energy = en_part(trial_particle, xi, yi, zi, xyz) - en_old;
    if (random->uniform() < exp(-beta*(delta_energy))) {
      energy += delta_energy;
      ++num_accepted;
      xyz[trial_particle][0] = xi;
      xyz[trial_particle][1] = yi;
      xyz[trial_particle][2] = zi;
    //} else {
      //xyz[trial_particle] = xyz_old;
    }
    energy_sum += energy;
    if (verbose) std::cout << "energy: " << energy << std::endl;
    if (trial % 10000 == 0) std::cout << "energy: " << energy << std::endl;
  }
  std::cout << "av en " << energy_sum / static_cast<double>(num_trials) << std::endl;
  std::cout << "acceptance " << num_accepted / static_cast<double>(num_trials) << std::endl;
}

// Conduct a parallel prefetch Metropolis Monte Carlo simulation.
void parallel(RandomBenchmark * random) {
  std::vector<std::vector<double> > xyz = read_xyz();
  std::uniform_real_distribution<double> ran_uniform;
  int proc_id, num_threads;
  #pragma omp parallel private(proc_id)
  {
    proc_id = omp_get_thread_num();
    if (proc_id == 0) {
      num_threads = static_cast<int>(omp_get_num_threads());
    }
  }
  double energy = en_all(xyz);
  double energy_sum = 0;
  int num_accepted = 0;
  int trial = 0;
  std::vector<bool> will_accept(num_threads);
  std::vector<double> delta_en(num_threads);
  std::vector<int> ppart(num_threads);
  std::vector<std::vector<double> > pxyz(num_threads, std::vector<double>(dimen));
  std::vector<RandomBenchmark> rans(num_threads);
  int first_accepted;
  #pragma omp parallel private(proc_id)
  {
    proc_id = omp_get_thread_num();
    rans[proc_id] = *random;
    if (proc_id == 0) {
      for (int i = 0; i < num_threads - 1; ++i) {
        for (int j = i + 1; j < num_threads; ++j) {
          rans[i].uniform();
        }
      }
    }
      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      rans[proc_id].seed();
      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

    #pragma omp barrier

    while (trial < num_trials) {

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      ppart[proc_id] = static_cast<int>(rans[proc_id].uniform()*num_particles);
      if (verbose) std::cout << "p " << ppart[proc_id] << std::endl;
      std::vector<double> * pixyz = &pxyz[proc_id];
      *pixyz = xyz[ppart[proc_id]];
      const double en_old = en_part(ppart[proc_id],
        (*pixyz)[0],
        (*pixyz)[1],
        (*pixyz)[2],
        xyz);
//      std::vector<double> xyz_old = xyz[ppart[proc_id]];
    double ran = rans[proc_id].uniform();
    if (verbose) std::cout << "ran " << ran <<std::endl;
    (*pixyz)[0] += max_move*(ran - 0.5);
      //(*pixyz)[0] += max_move*(rans[proc_id].uniform() - 0.5);
      (*pixyz)[1] += max_move*(rans[proc_id].uniform() - 0.5);
      (*pixyz)[2] += max_move*(rans[proc_id].uniform() - 0.5);
      if (verbose) std::cout << "new pos: " << (*pixyz)[0] << " " << (*pixyz)[1] << " " << (*pixyz)[2] << std::endl;
      delta_en[proc_id] = en_part(ppart[proc_id],
        (*pixyz)[0],
        (*pixyz)[1],
        (*pixyz)[2],
        xyz) - en_old;
      if (rans[proc_id].uniform() < exp(-beta*(delta_en[proc_id]))) {
        will_accept[proc_id] = true;
      } else {
        will_accept[proc_id] = false;
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #pragma omp barrier

      // determine first accepted and update Markov chain.
      #pragma omp critical
      {
        if (proc_id == 0) {
          first_accepted = num_threads;
          for (int i = 0; i < num_threads; ++i) {
            if (first_accepted == num_threads) {
              // haven't found accepted yet
              if (will_accept[i]) {
                first_accepted = i;
                energy += delta_en[i];
                xyz[ppart[i]] = pxyz[i];
                ++num_accepted;
                if (verbose) std::cout << "first_accepted " << first_accepted << std::endl;
              }
              energy_sum += energy;
              ++trial;
              if (verbose) std::cout << "energy: " << energy << std::endl;
              // std::cout << "trial " << trial << " energy " << energy << std::endl;
            }
          }
        }
      }
      #pragma omp barrier
    }
  }
  std::cout << "num trials " << trial << std::endl;
  std::cout << "av en " << energy_sum / static_cast<double>(trial) << std::endl;
  std::cout << "acceptance " << num_accepted / static_cast<double>(trial) << std::endl;
}

//int main() {
//  srand(time(NULL));  // seed by time
//  RandomBenchmark random;
//  // intialize and seed random number generator.
//  // std::cout << xyz[0][0] << std::endl;
//  // std::cout << xyz[182][2] << std::endl;
//  std::cout << en_all() << std::endl;
//  // std::cout << ran_uniform(generator) << std::endl;
//  // std::cout << ran_uniform(generator) << std::endl;
//  // std::cout << ran_uniform(generator) << std::endl;
//  // serial(&random);
//  parallel(&random);
//}

