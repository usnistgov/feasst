/**
This is an example usage of the C++ interface of FEASST to implement a
parallel prototype of the a single-component square well fluid.
In this parallel algorithm, a neighbor list is computed in parallel,
and then each trial is accepted/rejected in serial, with updates
to the neighbor list as needed.

The algorithm proceeds as follows:
 - For each thread, t (of n_t threads), initialize a copy
 - For each batch
   - For each thread in parallel
     - Attempt perturbation of Configuration and compute new map
   - Wait for all threads
   - For each thread, t, in serial
     - Decide to accept / reject.
     - If accepted
       - update configuration in all threads (synchronize) (in parallel, then wait)
       - update old neighbor list in all threads (including t) (in parallel, then wait)
       - update new neighbor list in remaining threads, t > n_t (in parallel, then wait)
 - End for each batch

The benefit of this parallel algorithm over Prefetch is that no trials are
"deleted" or "ignored" after the first accepted trial.
However, the overhead cost is higher than in Prefetch due to the neighbor list
updates required.
It turns out, both of these effects roughly cancel each other to yield roughly
the same efficiency as Prefetch.

Usage:
  mkdir build; cd &_
  cmake ..
  make
  ./main
*/

#include <iostream>
#include <unordered_set>
#include <vector>
#include "omp.h"
#include "feasst.h"

using namespace feasst;

static feasst::ArgumentParse args(
  "Canonical ensemble Prefetch Metropolis Monte Carlo simulation of Lennard Jones.", {
  {"--seed", "random number generator seed", "1234567"},
  {"--length", "cubic periodic boundary length", "9"},
  {"--disp", "trial displacement", "0.22"},
  //{"--disp", "trial displacement", "0.135"},
  {"--num", "number of particles", "500"},
  {"--data", "FEASST particle file", "/feasst/particle/lj.txt"},
  {"--beta", "inverse temperature", "1.2"},
  {"--trials", "number of Monte Carlo trials", "1e6"}});

const double cutoffsq = std::pow(1.5, 2);
const int dimen = 3;
//const int check = 1;
//const bool verbose = true;
//const int check = 1e6;
const int check = 1e4;
const bool verbose = false;

// Function to check if a vector contains any duplicate elements
template <typename T>
bool hasDuplicate(const std::vector<T>& vec) {
  std::unordered_set<T> seen;
  for (const auto& item : vec) {
    // If item already exists in the set, it's a duplicate
    if (!seen.insert(item).second) {
        return true;
    }
  }
  return false;
}

class Data {
 public:
  Data(const int num, const int dimen, const double length) : length_(length) {
    resize(num, dimen, &xyz);
    neighbor.resize(num);
    xold.resize(dimen);
  }
  std::vector<std::vector<double> > xyz;
  std::vector<std::vector<int> > neighbor;
  std::vector<int> new_neigh;
  std::vector<double> xold;
  int part;
  double entot;
  const double length_;
  bool overlap = false;
};

double energy(const int i, const int j, const std::vector<std::vector<double> >& xyz, const double length, double * r2, double * dx) {
  *r2 = 0;
  for (int dim = 0; dim < dimen; ++dim) {
    *dx = xyz[i][dim] - xyz[j][dim];
    if (*dx < -0.5*length) {
      *dx += length;
    } else if (*dx > 0.5*length) {
      *dx -= length;
    }
    *r2 += (*dx)*(*dx);
  }
  if (*r2 > cutoffsq) {
    return 0.;
  } else if (*r2 < 0.95) {
    return 1e250;
  } else {
    return -1;
  }
}

double energy(Data * data) {
  const int ipart = data->part;
  const double length = data->length_;
  const std::vector<std::vector<double> >& xyz = data->xyz;
  auto nwn = &(data->new_neigh);
  data->overlap = false;
  nwn->clear();
  double r2, dx, entot = 0;
  for (int jpart = 0; jpart < static_cast<int>(xyz.size()); ++jpart) {
    if (jpart != ipart) {
      const double en = energy(ipart, jpart, xyz, length, &r2, &dx);
      entot += en;
      if (en < -0.1) {
        nwn->push_back(jpart);
      } else if (en > 1) {
        data->overlap = true;
      }
    }
  }
  return entot;
}

double energy_all(Data * data) {
  double r2, dx, entot = 0.;
  const double length = data->length_;
  const int num = static_cast<int>(data->xyz.size());
  auto nb = &(data->neighbor);
  nb->clear();
  nb->resize(num);
  for (int i = 0; i < num - 1; ++i) {
    for (int j = i + 1; j < num; ++j) {
      const double en = energy(i, j, data->xyz, length, &r2, &dx);
      entot += en;
      if (en < -0.1) {
        (*nb)[i].push_back(j);
        (*nb)[j].push_back(i);
      }
    }
  }
  return entot;
}

void update_neighbor(const int part, const std::vector<int>& new_neigh, std::vector<std::vector<int> > * neighbor) {
  auto neigh = &(*neighbor)[part];
  // remove old neighbors
  for (int jpart : (*neigh)) {
    int findex = -1;
    ASSERT(find_in_list(part, (*neighbor)[jpart], &findex), "err");
    (*neighbor)[jpart].erase((*neighbor)[jpart].begin() + findex);
  }
  // add new neighbors
  for (int jpart : new_neigh) {
    auto nn = &(*neighbor)[jpart];
    nn->push_back(part);
    std::sort(nn->begin(), nn->end());
  }
  (*neigh) = new_neigh;
}

int main(int argc, char ** argv) {
  std::cout << "FEASST version " << feasst::FEASST_VERSION << std::endl
            << "# " << args.parse(argc, argv) << std::endl;
  const double length = args.get_double("--length");
  const int num = args.get_int("--num");
  const int batches = args.get_int("--trials");
  const double beta = args.get_double("--beta");
  const double disp = args.get_double("--disp");
  int num_threads = -1;
  #pragma omp parallel
  {
    num_threads = static_cast<int>(omp_get_num_threads());
  }
  std::vector<Data*> datum(num_threads);
  #pragma omp parallel
  {
    int t1 = omp_get_thread_num();
    if (verbose) {
      INFO("num:" << num);
      INFO("dimen:" << dimen);
    }
    Data data(num, dimen, length);
    datum[t1] = &data;
    double r2, dx;

    // init xyz
    std::ifstream file("../init.xyz");
    for (int i = 0; i < num; ++i) {
      for (int dim = 0; dim < dimen; ++dim) {
        double x;
        file >> x;
        datum[t1]->xyz[i][dim] = x;
      }
    }
    if (verbose) INFO("dat00 " << datum[t1]->xyz[0][0]);

    // init potential
    datum[t1]->entot = energy_all(datum[t1]);
    if (verbose) INFO("entot:" << datum[t1]->entot);
    ASSERT(datum[t1]->entot <= 0., datum[t1]->entot);

    Accumulator av_en;
    Accumulator accepted;
    std::unique_ptr<RandomMT19937> ran = std::make_unique<RandomMT19937>(argtype({{"seed", str(str_to_int(args.get("--seed"))+1234567*t1)}}));
      for (int batch = 0; batch < batches; ++batch) {
        #pragma omp barrier
//        if (t1 == 0) { for (int t2 = 0; t2 < num_threads; ++t2) {
        const int t2 = t1; {{
          if (verbose) INFO("**********BEGIN BATCH " << batch << " thread " << t2 << "**********");
          datum[t2]->part = ran->uniform(0, num - 1);
          if (verbose) INFO("part " << datum[t2]->part);
          datum[t2]->xold = datum[t2]->xyz[datum[t2]->part];
          for (int dim = 0; dim < dimen; ++dim) {
            datum[t2]->xyz[datum[t2]->part][dim] += (2.*ran->uniform() - 1.)*disp;
          }
          energy(datum[t2]);
        }}

        #pragma omp barrier
        if (t1 == 0) {
          for (int t2 = 0; t2 < num_threads; ++t2) {

//        #pragma omp for ordered schedule(static,1)
//        for (int t2 = 0; t2 < num_threads; ++t2) {
//          ASSERT(t1 == t2, "err");
//          #pragma omp ordered
//          {

            Data * dat = datum[t2];
            if (verbose) INFO("**********SYNC BATCH " << batch << " thread " << t2 << "**********");
            const int part = dat->part;
            double delta_en;
            if (dat->overlap) {
              delta_en = 1e250;
            } else {
              const double enold = -1.*static_cast<double>(dat->neighbor[part].size());
              if (verbose) INFO("enold " << enold);
              const double ennew = -1.*static_cast<double>(dat->new_neigh.size());
              if (verbose) INFO("ennew " << ennew);
              delta_en = ennew - enold;
            }
            if (verbose) INFO("delta_en " << delta_en);
            if (ran->uniform() < std::exp(-beta*delta_en)) {
              if (verbose) INFO("accepted");
              dat->entot += delta_en;
              accepted.accumulate(1);

              // update the neighbors of the current thread
              update_neighbor(part, dat->new_neigh, &(dat->neighbor));

              // update the position of all other threads
              for (int t3 = 0; t3 < num_threads; ++t3) {
                if (t3 != t2) {
                  datum[t3]->xyz[part] = dat->xyz[part];
                  if (part == datum[t3]->part) {
                    datum[t3]->xold = dat->xyz[part];
                  }
                }
              }

              // update map of all other threads
              for (int t3 = 0; t3 < num_threads; ++t3) {
                if (t3 != t2) {
                  // recompute energy if other particles have changed since computation
                  // or less efficiencly, always recompute
                  const double eij = energy(part, datum[t3]->part, datum[t3]->xyz, datum[t3]->length_, &r2, &dx);
                  if (verbose) INFO("eij " << eij);
                  datum[t3]->entot += delta_en;

                  update_neighbor(part, dat->new_neigh, &(datum[t3]->neighbor));
                  if (eij > 1) {
                    datum[t3]->overlap = true;
                  }
                  int findex = -1;
                  auto nn = &(datum[t3]->new_neigh);
                  const bool found = find_in_list(part, *nn, &findex);
                  if (eij < -0.1) {
                    // add to t3 nwn (if not there)
                    if (!found) {
                      nn->push_back(part);
                      std::sort(nn->begin(), nn->end());
                    }
                  } else {
                    // remove from t3 nwn (if there)
                    if (found) {
                      nn->erase(nn->begin() + findex);
                    }
                  }
                }
              }
            } else {
              if (verbose) INFO("rejected. " << feasst_str(dat->xold));
              dat->xyz[part] = dat->xold;
              accepted.accumulate(0);
            }
            if (verbose) INFO("entot: " << dat->entot);
            av_en.accumulate(dat->entot);

            // periodicly check that the xyz and map of each thread below t2 is equivalent
            if (batch % check == 0) {
            //if (batch % 1 == 0) {
              for (int t3 = 0; t3 < t2; ++t3) {
                if (verbose) INFO("checking xyz of t" << t3);
                for (int i = 0; i < num; ++i) {
                  for (int dim = 0; dim < dimen; ++dim) {
                    const double diff = dat->xyz[i][dim] - datum[t3]->xyz[i][dim];
                    if (std::abs(diff) > 1e-8) {
                      INFO(t2 << " " << i << " " << dim << " " << dat->xyz[i][dim]);
                      INFO(t3 << " " << i << " " << dim << " " << datum[t3]->xyz[i][dim]);
                      FATAL("i:" << i << " dim:" << dim << " diff:" << diff);
                    }
                  }
                }
              }
            }

            // periodically check for duplicate neighbor
            if (batch % check == 0) {
              for (int ip = 0; ip < num; ++ip) {
                if (hasDuplicate(dat->neighbor[ip])) {
                  FATAL(ip << " " << feasst_str(dat->neighbor[ip]));
                }
              }
            }

            // periodicly check that the recomputed neighbors agree
            // and that the total energy agrees with recomputed
            //if (batch % 1 == 0) {
            if (batch % check == 0) {
              std::vector<std::vector<int> > ncopy = dat->neighbor;
              //INFO(feasst_str(ncopy));
              energy_all(dat);
              ASSERT(is_equal(dat->neighbor, ncopy), "err");
              double ent = 0;
              for (int ip = 0; ip < num; ++ip) {
                ent -= static_cast<double>(dat->neighbor[ip].size())/2.;
              }
              ASSERT(std::abs(ent - dat->entot) < 1e-6, "ent:" << ent << " entot:" << dat->entot << " diff:" << ent-dat->entot);
            }
          }
        }
      }
      if (t1 == 0) {
        INFO("energy stat: " << av_en.str());
        INFO("acceptance stat: " << accepted.str());
        INFO(MAX_PRECISION << datum[t1]->entot);
    }
  }
}
