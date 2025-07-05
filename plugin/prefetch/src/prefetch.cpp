#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
#include <cmath>
#include <limits>
#include <random>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "threads/include/thread_omp.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/action.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "prefetch/include/prefetch.h"

// use this to make prefetch serial and simplify debugging
// #define DEBUG_SERIAL_MODE_5324634

namespace feasst {

Prefetch::Prefetch(argtype args) {
  activate_prefetch();
  trials_per_check_ = integer("trials_per_check", &args, 1e6);
  load_balance_ = boolean("load_balance", &args, false);
  ghost_ = boolean("ghost", &args, false);
  is_synchronize_ = boolean("synchronize", &args, false);
  #ifdef DEBUG_SERIAL_MODE_5324634
    WARN("DEBUG_SERIAL_MODE_5324634");
  #endif
  feasst_check_all_used(args);
}

void Prefetch::reset_trial_stats() {
  MonteCarlo::reset_trial_stats();
  for (Pool& pool : pool_) {
    pool.mc->reset_trial_stats();
  }
}

MonteCarlo * Prefetch::clone_(const int ithread) {
  if (ithread == 0) {
    return this;
  }
  return pool_[ithread].mc.get();
}

void Prefetch::create(std::vector<Pool> * pool) {
  // Initialize MC clones for each processor in pool_
  ASSERT(pool_.size() == 0, "pool is of size:" << pool_.size());
  if (ThreadOMP().is_enabled()) {
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #else
    {
    #endif
      num_threads_ = ThreadOMP().num();
    }
  } else {
    num_threads_ = 1;
  }
  pool_.resize(num_threads_);

  // set all trials for delayed finalization
  delay_finalize_();

  //clones_.resize(num_threads_);
  for (int thread = 1; thread < num_threads_; ++thread) {
    std::stringstream clone_ss;
    MonteCarlo::serialize(clone_ss);
    pool_[thread].mc = std::make_unique<MonteCarlo>(clone_ss);
  }

  // seed random number generators so that clones are not equal
  for (int i = 0; i < num_threads_ - 1; ++i) {
    clone_(i)->seed_random(rand());
  }

  // run some checks before attempting trials
  for (int thread = 0; thread < num_threads_; ++thread) {
    clone_(thread)->before_attempts_();
  }
}

void Prefetch::run_until_complete_(TrialFactory * trial_factory,
                                   Random * random) {
  DEBUG("here");
  if (!is_activated_) {
    MonteCarlo::run_until_complete_(trial_factory, random);
    return;
  }
  pool_.clear();
  attempt_(-1, trial_factory, random);
  write_checkpoint();
  write_to_file();
}

void Prefetch::attempt_(
    int num_trials,
    TrialFactory * trial_factory,
    Random * random) {
  DEBUG("Prefetch attempting");
  if (!is_activated_) {
    MonteCarlo::attempt_(num_trials, trial_factory, random);
    return;
  }
  if (num_trials < 10 && num_trials > 0) {
    WARN("inefficient use of prefetching with num_trials: " << num_trials
     << ". Consider using mc.attempt() instead of mc.run(MakeRun())");

  }

  // Require OPENMP; however, maintain ability to compile without.
  #ifndef _OPENMP
    FATAL("requires openmp");
  #endif // _OPENMP

  // If num_trials is -1, run based on criteria completion
  bool check_criteria_for_completion = false;
  if (num_trials == -1) {
    check_criteria_for_completion = true;
    num_trials = std::numeric_limits<int>::max();
  }

  if (pool_.size() == 0) {
    create(&pool_);
  }

  int itrial = 0;
  int first_thread_accepted;
  int proc_id = 0;
  #ifdef _OPENMP
  #pragma omp parallel private(proc_id)
  {
    proc_id = omp_get_thread_num();
  #endif // _OPENMP

    // num_trials may change to terminate loop
    bool complete = false;
    while (!complete) {
      #ifndef _OPENMP
      for (proc_id = 0; proc_id < num_threads_; ++proc_id) {
      #endif // _OPENMP

      //double previous_energy = criteria().current_energy();
      //DEBUG("previous_energy " << previous_energy);

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      if (proc_id == 0) {
        DEBUG("************************");
        DEBUG("* Begin Prefetch cycle *");
        DEBUG("************************");
        DEBUG("N " << clone_(proc_id)->configuration().num_particles());
      }

      // ordered list of index of trial to perform on each thread.
      if (proc_id == 0) {
        if (load_balance_) {
          // perform the same type of trial on each thread.
          const int index = trial_factory->random_index(random);
          for (int ithread = 0; ithread < num_threads_; ++ithread) {
            pool_[ithread].set_index(index);
          }
        } else {
          // randomly generator the type of trial for each thread.
          for (int ithread = 0; ithread < num_threads_; ++ithread) {
            pool_[ithread].set_index(trial_factory->random_index(random));
          }
        }
//        for (int ithread = 0; ithread < num_threads_; ++ithread) {
//          DEBUG("num attempts " << pool_[ithread].mc.trials().num_attempts() << " "
//                                << pool_[ithread].mc.trials().num_success() << " " << &pool_[ithread].mc);
//        }
//        DEBUG("num attempts main " << trials().num_attempts() << " "
//                                     << trials().num_success() << " " << this);
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      // set random number generators to store (zero storage for selecting trial?).
      Pool * pool = &pool_[proc_id];
      MonteCarlo * mc = clone_(proc_id);
      mc->load_cache_(true);

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      // Each processor attempts their trial in parallel,
      // without analyze modify or checkpoint.
      // Store new macrostate and acceptance prob
      pool->set_accepted(mc->attempt_trial(pool->index()));
      pool->set_auto_rejected(mc->trial(pool->index()).accept().reject());
      pool->set_endpoint(mc->trial(pool->index()).accept().endpoint());
      pool->set_ln_prob(mc->trial(pool->index()).accept().ln_metropolis_prob());
      DEBUG("proc id " << proc_id << " ln prob " << pool->ln_prob());
      DEBUG("critical proc_id " << proc_id << " " << pool->str());
      DEBUG("nump " << mc->system().configuration().num_particles());

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      // Determine first trial accepted (if any).
      if (proc_id == 0) {
        first_thread_accepted = num_threads_;
        for (int ithread = 0; ithread < num_threads_; ++ithread) {
          if (pool_[ithread].accepted() &&
              first_thread_accepted == num_threads_) {
            first_thread_accepted = ithread;
            DEBUG(MAX_PRECISION << ithread << " accepted, en: " << clone_(ithread)->criteria().current_energy());
          }
        }
        DEBUG("first thread " << first_thread_accepted);
      }

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      // any trial after accepted may contribute as a ghost
      if (proc_id == 0 && ghost_) {
        DEBUG("Update criteria with ghost trials (for TM) for each thread after first accepted");
        for (int ithread = first_thread_accepted + 1;
             ithread < num_threads_;
             ++ithread) {
          const Criteria& old_criteria = clone_(ithread)->criteria();
          for (int jthread = 0; jthread < num_threads_; ++jthread) {
            DEBUG("ghost " << ithread << " " << jthread);
            DEBUG("first accepted " << first_thread_accepted);
            clone_(jthread)->ghost_trial_(
              pool_[ithread].ln_prob(),
              old_criteria.state_old(),
              old_criteria.state_new(),
              pool_[ithread].endpoint());
          }
        }
      }

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      // revert trials after accepted trial.
      if (first_thread_accepted != num_threads_) {
        if (proc_id > first_thread_accepted) {
          DEBUG("reverting trial " << proc_id);
          mc->revert_(pool->index(), pool->accepted(), pool->endpoint(),
                      pool->auto_rejected(), pool->ln_prob());
        }
      }


      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      // testing decouple here HWH
      if (proc_id == 0) {
        DEBUG("for each thread up to the first accepted, "
           << "update other threads (incl. main) regarding failed attempt by "
           << "thread. update steppers");
        for (int ithread = 0; ithread < first_thread_accepted; ++ithread) {
          DEBUG("imitate failed " << ithread);
          // loop through other threads to update
          const Criteria& old_criteria = clone_(ithread)->criteria();
          DEBUG("nump " << clone_(ithread)->configuration().num_particles());
          for (int jthread = 0; jthread < num_threads_; ++jthread) {
            DEBUG("j " << jthread);
            if (ithread != jthread) {
              DEBUG(" index " << pool_[ithread].index() <<
                   " lnp " << pool_[ithread].ln_prob() <<
                   " old " << old_criteria.state_old() <<
                   " new " << old_criteria.state_new());
              clone_(jthread)->imitate_trial_rejection_(
                pool_[ithread].index(),
                pool_[ithread].ln_prob(),
                pool_[ithread].endpoint(),
                pool_[ithread].auto_rejected(),
                old_criteria.state_old(),
                old_criteria.state_new()
              );
            } else {
              DEBUG("Update TM on rejection");
              DEBUG(clone_(jthread)->get_criteria()->class_name());
              clone_(jthread)->get_criteria()->imitate_trial_rejection_(
                pool_[ithread].ln_prob(),
                old_criteria.state_old(),
                old_criteria.state_new(),
                pool_[ithread].endpoint()
              );
              DEBUG("here");
            }
            if (jthread == 0) {
              DEBUG("analyzing");
              after_trial_analyze_();
              DEBUG("here");
            }
          }
//          DEBUG("num attempts " << pool->mc->trials().num_attempts() << " "
//                                << pool->mc->trials().num_success() << " " << &pool->mc);
        }
        DEBUG("num attempts main " << trials().num_attempts() << " "
                                     << trials().num_success() << " " << this);
      }

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      if (first_thread_accepted < num_threads_) {
        DEBUG("Replicate first accepted trial in all other threads in proc_id " << proc_id);
        Pool * accepted_pool = &pool_[first_thread_accepted];
        if (proc_id != first_thread_accepted) {
          // load/unload system energies and random numbers
          mc->unload_cache_(*clone_(first_thread_accepted));
          mc->attempt_trial(accepted_pool->index());
        }
        mc->finalize_(accepted_pool->index());
      } else {
        DEBUG("all rejected, en: " << criteria().current_energy());
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      if (first_thread_accepted < num_threads_) {
        if (is_synchronize_) {
          DEBUG("synchronize other threads with first accepted thread " << proc_id);
          Pool * accepted_pool = &pool_[first_thread_accepted];
          if (proc_id != first_thread_accepted) {
            DEBUG(proc_id);
            DEBUG(first_thread_accepted);
            DEBUG(accepted_pool->index());
            const MonteCarlo& cln = *clone_(first_thread_accepted);
            DEBUG(cln.trial(accepted_pool->index()).accept().perturbed(0).str());
            mc->synchronize_(cln,
              cln.trial(accepted_pool->index()).accept().perturbed());
          }
        }
        if (proc_id == 0) {
          after_trial_analyze_();
        }
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      // disable cache
      mc->load_cache_(false);

      // update last trial for tuning
      int last_thread = first_thread_accepted;
      if (last_thread >= num_threads_) {
        last_thread = num_threads_ - 1;
      }
      if (proc_id != last_thread) {
        const MonteCarlo& cln = *clone_(last_thread);
        mc->get_trial_factory()->set_last_index(cln.trials().last_index());
        mc->get_criteria()->set_was_accepted(cln.criteria().was_accepted());
      }

      // perform after trial on all clones/main after multiple trials performed
      // do this in serial so that files are not written to by multiple threads
      // simultaneously
      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #pragma omp critical
      {
      #else
      {
      #endif
      for (int im = 0;
           im < std::min(num_threads_, first_thread_accepted + 1);
           ++im) {
        // DEBUG("im " << im << " first " << first_thread_accepted);
        mc->after_trial_modify_();
        mc->after_trial_checkpoint_();
      }
      }

      DEBUG("update itrial");
      //if (proc_id <= first_thread_accepted) {
      //  ++itrial;
      //  ++trials_since_check_;
      if (proc_id == 0) {
        int increment = num_threads_;
        if (first_thread_accepted < num_threads_) {
          increment = first_thread_accepted + 1;
        }
        itrial += increment;
        trials_since_check_ += increment;
      }

//      #ifdef DEBUG_SERIAL_MODE_5324634
//      }
//      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      DEBUG("periodically check that all threads are equal");
      if (trials_since_check_ >= trials_per_check_ && proc_id > 0) {
        trials_since_check_ = 0;
        const double tolerance = 1e-8;
        const MonteCarlo& mcc = *clone_(proc_id);
        for (int conf = 0; conf < system().num_configurations(); ++conf) {
          DEBUG("conf " << conf);
          const double energy = criteria().current_energy(conf);
          DEBUG("check that the current energy of all threads and main are the same: " << energy);
          const double diff = mcc.criteria().current_energy(conf) - energy;
          ASSERT(std::abs(diff) <= tolerance, "diff: " << diff);
          DEBUG("check that conf:" << conf << " on proc_id: " << proc_id << " is equal to the first proc");
          ASSERT(system().configuration(conf).is_equal(mcc.system().configuration(conf), tolerance), "configs not equal thread" << proc_id);
        }
        ASSERT(trials().is_equal(mcc.trials()), "trials not equal thread" << proc_id);
        ASSERT(criteria().is_equal(mcc.criteria(), tolerance), "criteria not equal: " << proc_id);

          // periodically equate all clones to first thread to prevent drift
  //        std::stringstream clone_ss;
  //        MonteCarlo::serialize(clone_ss);
  //        pool_[proc_id].mc = MonteCarlo(clone_ss);
  //        clone_(proc_id)->seed_random(rand());  // HWH not thread safe // HWH consider setting RNG to previous seeds manually
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      DEBUG("itrial: " << itrial);
      if (check_criteria_for_completion) {
        if (criteria().is_complete()) {
          complete = true;
        }
      } else {
        if (itrial >= num_trials) {
          complete = true;
        }
      }
      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP
    }
  }
}

void Prefetch::run(std::shared_ptr<Action> action) {
  DEBUG("is_activated_ " << is_activated_);
  DEBUG("action class name: " << action->class_name());
  if (is_activated_ && static_cast<int>(pool_.size()) > 0) {
    DEBUG("here");
    if (action->class_name() == "Run") {
      action->run(this);
    } else {
      DEBUG("here " << pool_.size());
      for (int thread = 0; thread < num_threads_; ++thread) {
        MonteCarlo * mc = clone_(thread);
        mc->MonteCarlo::run(action);
      }
    }
  } else {
    MonteCarlo::run(action);
  }
}

void Prefetch::run_num_trials(int64_t num_trials) {
  if (num_trials < 0) {
    return;
  }
  if (num_trials < 10) {
    FATAL("num_trials:" << num_trials << " is inefficient.");
  }
  pool_.clear();
  while (num_trials > 0) {
    attempt(num_trials);
  }
}

void Prefetch::run_until_num_particles(const int num_particles,
                                       const std::string& particle_type,
                                       const int configuration_index) {
  activate_prefetch(false);
  MonteCarlo::run_until_num_particles(num_particles, particle_type,
                                      configuration_index);
  activate_prefetch(true);
}

void Prefetch::run_until_file_exists(const std::string& file_name,
    const int trials_per_file_check) {
  if (!file_name.empty()) {
    WARN("run_until_file is not implemented efficiently with Prefetch.");
  }
  MonteCarlo::run_until_file_exists(file_name, trials_per_file_check);
}

void Prefetch::serialize(std::ostream& ostr) const {
  MonteCarlo::serialize(ostr);
  feasst_serialize_version(5686, ostr);
  feasst_serialize(is_activated_, ostr);
  feasst_serialize(trials_per_check_, ostr);
  feasst_serialize(trials_since_check_, ostr);
  feasst_serialize(load_balance_, ostr);
  feasst_serialize(is_synchronize_, ostr);
  feasst_serialize(ghost_, ostr);
}

Prefetch::Prefetch(std::istream& istr) : MonteCarlo(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5686, "version: " << version);
  feasst_deserialize(&is_activated_, istr);
  feasst_deserialize(&trials_per_check_, istr);
  feasst_deserialize(&trials_since_check_, istr);
  feasst_deserialize(&load_balance_, istr);
  feasst_deserialize(&is_synchronize_, istr);
  feasst_deserialize(&ghost_, istr);
}

const std::string Pool::str() const {
  std::stringstream ss;
  ss << index_ << " " << ln_prob_ << " " << accepted_;
  return ss.str();
}

void Prefetch::run_until_complete() {
  DEBUG("here");
  run_until_complete_(get_trial_factory(), get_random());
}

}  // namespace feasst
