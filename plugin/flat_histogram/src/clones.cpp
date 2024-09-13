#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
//#ifdef _WIN32
//  #include <Windows.h>  // sleep
//#else
//  #include <unistd.h>  // sleep
//#endif
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <fstream>
#include "utils/include/custom_exception.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/checkpoint.h"
#include "math/include/constants.h"
#include "math/include/histogram.h"
#include "math/include/utils_math.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/clones.h"

namespace feasst {

Clones::~Clones() {}

//void Clones::create(const MonteCarlo& mc,
//                    const Window& window) {
//  std::stringstream ss_fh, ss_macro;
//  ASSERT(mc.criteria().class_name() == "FlatHistogram",
//    "assumes flat histogram criteria: " << mc.criteria().class_name());
//  mc.criteria().serialize(ss_fh);
//  FlatHistogram fh(ss_fh);
//  fh.macrostate().serialize(ss_macro);
//  auto macrostate = std::make_shared<Macrostate>(ss_macro);
//  for (std::vector<int> bounds : window.boundaries()) {
//    auto clone = std::make_shared<MonteCarlo>(deep_copy(mc));
//    auto macro2 = deep_copy_derived(macrostate);
//    macro2->set(Histogram({{"width"
//    auto macrostate =
//    clone->set(MakeFlatHistogram(
//  }
//}

const MonteCarlo& Clones::clone(const int index) const {
  ASSERT(index < num(), "index: " << index << " >= num: " << num());
  return const_cast<MonteCarlo&>(*clones_[index]);
}

MonteCarlo * Clones::get_clone(const int index) {
  ASSERT(index < num(), "index: " << index << " >= num: " << num());
  return clones_[index].get();
}

void Clones::initialize(const int upper_index, argtype args) {
  const int attempt_batch = integer("attempt_batch", &args, 1);
  const int max_batch = integer("max_batch", &args, -1);
  feasst_check_all_used(args);
  Acceptance empty;
  MonteCarlo * lower = clones_[upper_index - 1].get();
  MonteCarlo * upper = clones_[upper_index].get();
  std::unique_ptr<FlatHistogram> fh_lower = flat_histogram(upper_index - 1);
  std::unique_ptr<FlatHistogram> fh_upper = flat_histogram(upper_index);
  const double macro_upper_min = fh_upper->macrostate().value(0);
  const double macro_lower_max =
    fh_lower->macrostate().histogram().center_of_last_bin();

  if (fh_upper->macrostate().is_allowed(upper->system(),
                                        upper->criteria(), empty)) {
    DEBUG("already initialized");
    return;
  }

  DEBUG("find configuration in lower which overlaps with upper");
  bool overlap = false;
  int batch = 0;
  while (!overlap) {
    lower->attempt(attempt_batch);

    const double macro_lower = fh_lower->macrostate().value(
      lower->system(), lower->criteria(), empty);
    DEBUG("macro_lower: " << macro_lower);
    if (is_in_interval(macro_lower, macro_upper_min, macro_lower_max)) {
      overlap = true;
    }

    ++batch;
    ASSERT(max_batch == -1 || batch < max_batch,
      "reached maximum batch: " << batch);
  }

  DEBUG("send configuration to upper");
  //upper->set(deep_copy(lower->system()));
  upper->get_system()->get_configuration()->copy_particles(
    lower->configuration(),
    true);
  //std::cout << "Info 1 num " << upper->configuration().num_particles() << std::endl;
  DEBUG("num " << upper->configuration().num_particles());
  upper->initialize_criteria();
}

void Clones::initialize(argtype args) {
  for (int upper_index = 1; upper_index < num(); ++upper_index) {
    initialize(upper_index, args);
  }
}

void Clones::run_until_complete(argtype args) {
#ifdef _OPENMP
  run_until_complete_omp_(args);
#else // _OPENMP
  run_until_complete_serial_();
#endif // _OPENMP
}

void Clones::run_until_complete_serial_() {
  DEBUG("run_until_complete_serial_");
  for (auto clone : clones_) {
    clone->run_until_complete();
  }
}

bool are_all_complete(const std::vector<bool>& is_complete) {
  bool complete = true;
  for (const bool comp : is_complete) {
    if (!comp) complete = false;
  }
  return complete;
}

void Clones::run_until_complete_omp_(argtype run_args,
                                     const bool init,
                                     argtype init_args) {
  DEBUG("run_until_complete_omp_");
#ifdef _OPENMP
  const int omp_batch = integer("omp_batch", &run_args, 1e6);
  std::string ln_prob_file;
  if (used("ln_prob_file", run_args)) {
    ln_prob_file = str("ln_prob_file", &run_args);
  }
  feasst_check_all_used(run_args);
  std::vector<bool> is_complete(num(), false);
  std::vector<bool> is_initialized(num(), false);
  is_initialized[0] = true;
  #pragma omp parallel
  {
    const int thread = omp_get_thread_num();
    const int num_thread = static_cast<int>(omp_get_num_threads());
    DEBUG("thread " << thread << " of " << num_thread);

    bool terminated = false;
    try {
      ASSERT(num() <= num_thread, "more clones: " << num() << " than OMP threads:"
        << num_thread << ". Use \"export OMP_NUM_THREADS=\" to set OMP threads.");
      if (thread < num()) {
        MonteCarlo * clone = clones_[thread].get();
        if (init) {
          // wait until thread is initialized
          while (!is_initialized[thread]) {
            std::this_thread::sleep_for (std::chrono::seconds(1));
          }
          //while (!is_initialized[thread]) sleep(0.01);
          DEBUG("thread " << thread << " is initialized");

          // initialize next thread, if applicable
          if (thread < num() - 1) initialize(thread + 1, init_args);
          is_initialized[thread + 1] = true;
        }

        // continue running while waiting for all threads to complete
        if (clone->criteria().is_complete()) is_complete[thread] = true;
        while (!are_all_complete(is_complete)) {
          clone->attempt(omp_batch);
          if (clone->criteria().is_complete()) is_complete[thread] = true;
          if (thread == 0) {
            if (!ln_prob_file.empty()) {
              std::ofstream file;
              file.open(ln_prob_file);
              for (const double value : ln_prob().values()) {
                file << value << std::endl;
              }
              file.close();
            }
          }
        }
      }
    } catch(const feasst::CustomException& e) {
      WARN(e.what());
      terminated = true;
    }
    DEBUG("terminated: " << terminated);

    if (thread == 0 && checkpoint_) checkpoint_->write(*this);
    if (thread < num()) clones_[thread]->write_checkpoint();

    #pragma omp barrier
    if (terminated && thread == 0) {
      FATAL("Clones::run_until_complete_omp was terminated.");
    }
  }

#else // _OPENMP
FATAL("Not complied with OMP");
#endif // _OPENMP
}

void Clones::initialize_and_run_until_complete(argtype run_args,
                                               argtype init_args) {
#ifdef _OPENMP
  run_until_complete_omp_(run_args, true, init_args);
#else // _OPENMP
  initialize(init_args);
  run_until_complete_serial_();
#endif // _OPENMP
}

std::unique_ptr<FlatHistogram> Clones::flat_histogram(const int index) const {
  std::stringstream ss;
  clone(index).criteria().serialize(ss);
  return std::make_unique<FlatHistogram>(ss);
}

LnProbability Clones::ln_prob(Histogram * macrostates,
    std::vector<double> * multistate_data,
    const std::string analyze_name,
    const AnalyzeData& get) const {
  std::vector<double> ln_prob;
  std::vector<double> edges;
  double shift = 0.;
  int starting_lower_bin = 0;
  for (int fh_index = 0; fh_index < num() - 1; ++fh_index) {
    DEBUG("fh_index " << fh_index);
    std::unique_ptr<FlatHistogram> fh_lower = flat_histogram(fh_index);
    std::unique_ptr<FlatHistogram> fh_upper = flat_histogram(fh_index + 1);
    const double macro_upper_min = fh_upper->macrostate().value(0);

    // Optionally, extract multistate_data
    std::vector<double> lower_data, upper_data;
    if (multistate_data) {
      SeekAnalyze seek;
      lower_data = seek.multistate_data(analyze_name, clone(fh_index), get);
      upper_data = seek.multistate_data(analyze_name, clone(fh_index + 1), get);
    }

    int upper_index = 0;
    std::vector<double> overlap_upper;
    std::vector<double> overlap_lower;
    DEBUG("starting_lower_bin " << starting_lower_bin);
    DEBUG("lower size " << fh_lower->bias().ln_prob().size());
    for (int bin = starting_lower_bin;
         bin < fh_lower->bias().ln_prob().size();
         ++bin) {
      DEBUG("bin " << bin);
      const double macro_lower = fh_lower->macrostate().value(bin);
      const double ln_prob_lower = fh_lower->bias().ln_prob().value(bin);
      if (std::abs(macro_lower - macro_upper_min) > NEAR_ZERO &&
          upper_index == 0) {
        DEBUG("macro_lower " << macro_lower);
        DEBUG("macro_upper_min " << macro_upper_min);
        ln_prob.push_back(ln_prob_lower + shift);
        if (multistate_data) multistate_data->push_back(lower_data[bin]);
      } else {
        DEBUG("upper index " << upper_index);
        DEBUG("stitch " << bin);
        overlap_upper.push_back(fh_upper->bias().ln_prob().value(upper_index));
        overlap_lower.push_back(ln_prob_lower);
        DEBUG("upper " << fh_upper->bias().ln_prob().value(upper_index));
        DEBUG("lower " << ln_prob_lower);
        if (multistate_data) {
          DEBUG("lower data " << lower_data[bin]);
          DEBUG("upper data " << upper_data[upper_index]);
          multistate_data->push_back(
            0.5*(lower_data[bin] + upper_data[upper_index]));
        }
        ++upper_index;
      }
      if (macrostates) {
        const double lower = fh_lower->macrostate().histogram().edges()[bin];
        edges.push_back(lower);
      }
    }
    ASSERT(upper_index > 0, "No overlap when upper_index: " << upper_index);

    // average the ln_probs of the overlapping region
    Accumulator ln_prob_shift;
    for (int ovl = 0; ovl < static_cast<int>(overlap_upper.size()); ++ovl) {
      ln_prob_shift.accumulate(overlap_lower[ovl] - overlap_upper[ovl]);
    }
    DEBUG("av shift: " << ln_prob_shift.average());
    for (int ovl = 0; ovl < static_cast<int>(overlap_upper.size()); ++ovl) {
      ln_prob.push_back(0.5*(overlap_lower[ovl] + overlap_upper[ovl] + ln_prob_shift.average()) + shift);
    }
    shift += ln_prob_shift.average();
    DEBUG("total shift " << shift);
    starting_lower_bin = static_cast<int>(overlap_lower.size());
  }

  // now add the non-overlapping part of the last clone
  std::unique_ptr<FlatHistogram> fh = flat_histogram(num() - 1);
  std::vector<double> data;
  if (multistate_data) {
    data = SeekAnalyze().multistate_data(analyze_name, clone(num() - 1), get);
  }
  for (int bin = starting_lower_bin; bin < fh->bias().ln_prob().size(); ++bin) {
    ln_prob.push_back(fh->bias().ln_prob().value(bin) + shift);
    if (macrostates) {
      const double macro = fh->macrostate().histogram().edges()[bin];
      edges.push_back(macro);
    }
    if (multistate_data) multistate_data->push_back(data[bin]);
  }
  if (macrostates) {
    edges.push_back(fh->macrostate().histogram().edges().back());
    *macrostates = Histogram();
    macrostates->set_edges(edges);
  }

  LnProbability lnpi(ln_prob);
  lnpi.normalize();
  return lnpi;
}

void Clones::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2845, ostr);
  feasst_serialize(clones_, ostr);
//  feasst_serialize(checkpoint_, ostr);
  feasst_serialize_endcap("Clones", ostr);
}

Clones::Clones(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2845, "version: " << version);
  // HWH for unknown reasons, this does not work
  //feasst_deserialize(&clones_, istr);
  int dim1;
  istr >> dim1;
  clones_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      clones_[index] = std::make_shared<MonteCarlo>(istr);
    }
  }
//  // HWH for unknown reasons, this function template does not work.
//  //feasst_deserialize(checkpoint_, istr);
//  { int existing;
//    istr >> existing;
//    if (existing != 0) {
//      checkpoint_ = std::make_shared<Checkpoint>(istr);
//    }
//  }
  feasst_deserialize_endcap("Clones", istr);
}

void Clones::set(std::shared_ptr<Checkpoint> checkpoint) {
  checkpoint_ = checkpoint;
}

std::shared_ptr<Clones> MakeClones(const std::vector<std::string> file_names) {
  auto clones = std::make_shared<Clones>();
  for (const std::string& file_name : file_names) {
    std::ifstream file(file_name);
    std::string line;
    std::getline(file, line);
    ASSERT(!line.empty(), "file: " << file_name << " is empty.");
    std::stringstream ss(line);
    clones->add(std::make_shared<MonteCarlo>(ss));
  }
  return clones;
}

std::shared_ptr<Clones> MakeClones(const std::string prepend,
  const int num,
  const int min,
  const std::string append) {
  std::vector<std::string> file_names;
  for (int index = min; index < num - min; ++index) {
    file_names.push_back(prepend + str(index) + append);
  }
  return MakeClones(file_names);
}

void Clones::set_num_iterations_to_complete(const int iterations) {
  for (int index = 0; index < num(); ++index) {
    clones_[index]->get_criteria()->set_num_iterations_to_complete(iterations);
  }
}

std::string Clones::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Clones Clones::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Clones(ss);
}

}  // namespace feasst
