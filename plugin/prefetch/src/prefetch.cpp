
// HWH consider making a feasst omp plugin that could interface other parallelization methods
#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
#include <cmath>
#include <limits>
#include "utils/include/serialize.h"
#include "prefetch/include/prefetch.h"

// // use this to turn prefetch serial to simply debugging
#define DEBUG_SERIAL_MODE_5324634

namespace feasst {

MonteCarlo * Prefetch::clone_(const int ithread) {
  if (ithread == 0) {
    return this;
  }
  return &pool_[ithread].mc;
}

void Prefetch::run_until_complete_(TrialFactory * trial_factory,
                                   Random * random) {
  if (!is_activated_) {
    MonteCarlo::run_until_complete_(trial_factory, random);
    return;
  }
  attempt_(-1, trial_factory, random);
}

void Prefetch::attempt_(
    int num_trials,
    TrialFactory * trial_factory,
    Random * random) {
  if (!is_activated_) {
    MonteCarlo::attempt_(num_trials, trial_factory, random);
    return;
  }

  // If num_trials is -1, run based on criteria completion
  bool check_criteria_for_completion = false;
  if (num_trials == -1) {
    check_criteria_for_completion = true;
    num_trials = std::numeric_limits<int>::max();
  }

  #ifndef _OPENMP
    ERROR("requires openmp");
  #endif // _OPENMP

  // initialize MC clones for each processor
  // and pool size.
  // alternatively, input the clones by using a functional
  // approach to creating mc objects
  int proc_id;
  int first_thread_accepted;
  double sum_en = 0;
  if (pool_.size() == 0) {
    #ifdef _OPENMP
      #pragma omp parallel private(proc_id)
      {
        proc_id = omp_get_thread_num();
        DEBUG("hello from " << proc_id);
        if (proc_id == 0) {
          num_threads_ = static_cast<int>(omp_get_num_threads());
        }
      }
    #else // _OPENMP
      num_threads_ = 1;
    #endif // _OPENMP
//    num_threads_ = 2;
    pool_.resize(num_threads_);

    // set all trials for delayed finalization
    delay_finalize();

    //clones_.resize(num_threads_);
    for (int thread = 1; thread < num_threads_; ++thread) {
      std::stringstream clone_ss;
      MonteCarlo::serialize(clone_ss);
      pool_[thread].mc = MonteCarlo(clone_ss);
    }

    // offset random number chains so that clones are not equal
    for (int i = 0; i < num_threads_ - 1; ++i) {
      for (int j = i + 1; j < num_threads_; ++j) {
        clone_(i)->offset_random();
      }
    }

    // run some checks before attempting trials
    for (int thread = 0; thread < num_threads_; ++thread) {
      clone_(thread)->before_attempts_();
    }
  }

  int itrial = 0;
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

      if (proc_id == 0) {
        DEBUG("************************");
        DEBUG("* Begin Prefetch cycle *");
        DEBUG("************************");
      }

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      // ordered list of randomly generated moves identified by index of trial.
      if (proc_id == 0) {
        if (load_balance_) {
          const int index = trial_factory->random_index(random);
          for (int ithread = 0; ithread < num_threads_; ++ithread) {
            pool_[ithread].set_index(index);
          }
        } else {
          for (int ithread = 0; ithread < num_threads_; ++ithread) {
            pool_[ithread].set_index(trial_factory->random_index(random));
          }
        }
        for (int ithread = 0; ithread < num_threads_; ++ithread) {
          DEBUG("num attempts " << pool_[ithread].mc.trials().num_attempts() << " "
                                << pool_[ithread].mc.trials().num_success() << " " << &pool_[ithread].mc);
        }
        DEBUG("num attempts master " << trials().num_attempts() << " "
                                     << trials().num_success() << " " << this);
      }

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      // set random number generators to store (zero storage for selecting trial?).
      Pool * pool = &pool_[proc_id];
      MonteCarlo * mc = clone_(proc_id);
      mc->load_cache(true);

      #ifdef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #endif

      // Each processor attempts their trial in parallel,
      // without analyze modify or checkpoint.
      // Store new macrostate and acceptance prob
      pool->set_accepted(mc->attempt_trial(pool->index()));
      pool->set_ln_prob(mc->trial(pool->index())->accept().ln_metropolis_prob());
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
            INFO(MAX_PRECISION << ithread << " accepted, en: " << clone_(ithread)->criteria()->current_energy());
          }
        }
        DEBUG("first thread " << first_thread_accepted);
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
          mc->revert(pool->index(), pool->accepted(), pool->ln_prob());
        }
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      #ifdef _OPENMP
      #pragma omp barrier
      #endif // _OPENMP

      if (proc_id == 0) {
        DEBUG("for each thread up to the first accepted, "
           << "update other threads (incl. main) regarding failed attempt by "
           << "thread. update steppers");
        for (int ithread = 0; ithread < first_thread_accepted; ++ithread) {
          DEBUG("imitate failed " << ithread);
          // loop through other threads to update
          const Criteria * old_criteria = clone_(ithread)->criteria();
          DEBUG("nump " << clone_(ithread)->configuration().num_particles());
          for (int jthread = 0; jthread < num_threads_; ++jthread) {
            if (ithread != jthread) {
              DEBUG("j " << jthread << " lnp " << pool_[ithread].ln_prob() <<
                   " old " << old_criteria->state_old() <<
                   " new " << old_criteria->state_new());
              clone_(jthread)->imitate_trial_rejection(
                pool_[ithread].index(),
                pool_[ithread].ln_prob(),
                old_criteria->state_old(),
                old_criteria->state_new()
              );
            }
            if (jthread == 0) {
              clone_(jthread)->after_trial_analyze_();
            }
            if (jthread == 0) {
              sum_en += clone_(jthread)->criteria()->current_energy();
            }
          }
          DEBUG("num attempts " << pool->mc.trials().num_attempts() << " "
                                << pool->mc.trials().num_success() << " " << &pool->mc);
        }
        DEBUG("num attempts master " << trials().num_attempts() << " "
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
        DEBUG("Replicate first accepted trial in all other threads");
        Pool * accepted_pool = &pool_[first_thread_accepted];
        if (proc_id != first_thread_accepted) {
          // load/unload system energies and random numbers
          mc->unload_cache(*clone_(first_thread_accepted));
          mc->attempt_trial(accepted_pool->index());
        }
        mc->finalize(accepted_pool->index());
        if (proc_id == 0) {
          after_trial_analyze_();
        }
        if (proc_id == 0) {
          sum_en += mc->criteria()->current_energy();
        }
      } else {
        INFO("all rejected, en: " << criteria()->current_energy());
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
      mc->load_cache(false);

      // perform after trial on all clones/master after multiple trials performed
      // do this in serial so that files are not written to by multiple threads
      // simultaneously
      #ifndef DEBUG_SERIAL_MODE_5324634
      #pragma omp critical
      {
      #else
      {
      #endif
      mc->after_trial_modify_();
      }

      DEBUG("update itrial");
      if (proc_id <= first_thread_accepted) {
        ++itrial;
        ++steps_since_check_;
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

      DEBUG("periodically check that all threads are equal");
      if (steps_since_check_ >= steps_per_check_ && proc_id > 0) {
        steps_since_check_ = 0;
        const double energy = criteria()->current_energy();
        DEBUG("check that the current energy of all threads and master are the same: " << energy);
        const double tolerance = 1e-8;
        const MonteCarlo& mcc = *clone_(proc_id);
        const double diff = mcc.criteria()->current_energy() - energy;
        ASSERT(std::fabs(diff) <= tolerance, "diff: " << diff);
        ASSERT(system().configuration().is_equal(mcc.system().configuration(), tolerance), "configs not equal thread" << proc_id);
        ASSERT(trials().is_equal(mcc.trials()), "trials not equal thread" << proc_id);
//        ASSERT(criteria()->is_equal(mcc.criteria()), "criteria not equal: " << proc_id);

        // periodically equate all clones to first thread to prevent drift
        std::stringstream clone_ss;
        MonteCarlo::serialize(clone_ss);
        pool_[proc_id].mc = MonteCarlo(clone_ss);
      }

      #ifdef DEBUG_SERIAL_MODE_5324634
      }
      #endif

      DEBUG("itrial: " << itrial);
      if (check_criteria_for_completion) {
        if (criteria()->is_complete()) {
          complete = true;
        }
      } else {
        if (itrial >= num_trials) {
          complete = true;
        }
      }
      if (proc_id == 0) {
        INFO("sum_en: " << sum_en << " n " << itrial << " av " << sum_en/itrial);
      }
    }
  }
}

void Prefetch::serialize(std::ostream& ostr) const {
  MonteCarlo::serialize(ostr);
  feasst_serialize_version(348, ostr);
  feasst_serialize(is_activated_, ostr);
  feasst_serialize(steps_per_check_, ostr);
  feasst_serialize(steps_since_check_, ostr);
  feasst_serialize(load_balance_, ostr);
}

Prefetch::Prefetch(std::istream& istr) : MonteCarlo(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 348, "version: " << version);
  feasst_deserialize(&is_activated_, istr);
  feasst_deserialize(&steps_per_check_, istr);
  feasst_deserialize(&steps_since_check_, istr);
  feasst_deserialize(&load_balance_, istr);
}

}  // namespace feasst
