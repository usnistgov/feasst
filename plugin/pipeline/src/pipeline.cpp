
#include <omp.h>
#include "pipeline/include/pipeline.h"

namespace feasst {

//void Pipeline::distribute_for_() {
//}

MonteCarlo * Pipeline::clone_(const int ithread, const int jthread) {
  if (ithread == jthread) {
    // update main
    return this;
  }
  return &clones_[jthread];
}

void Pipeline::attempt_(
    System * system,
    Criteria * criteria,
    TrialFactory * trial_factory,
    AnalyzeFactory * analyze_factory,
    ModifyFactory * modify_factory,
    Checkpoint * checkpoint,
    Random * random) {

  DEBUG("************************");
  DEBUG("* Begin Pipeline cycle *");
  DEBUG("************************");

  // initialize MC clones for each processor
  // and pool size.
  int proc_id;
  if (clones_.size() == 0) {
    #pragma omp parallel private(proc_id)
    {
      proc_id = omp_get_thread_num();
      INFO("hello from " << proc_id);
      if (proc_id == 0) {
        num_threads_ = static_cast<int>(omp_get_num_threads());
      }
    }
    clones_.resize(num_threads_);
    for (MonteCarlo& clone : clones_) {
      clone = *this;
    }
    trial_pool_.resize(num_threads_);
  }

  // Generate ordered list of randomly generated moves identified by index of
  // trial.
  INFO("num att " << trials().num_attempts());
  for (int ithread = 0; ithread < num_threads_; ++ithread) {
    trial_pool_[ithread].set_index(trial_factory->random_index(random));
    INFO("num attempts " << clones_[ithread].trials().num_attempts());
  }

  // Each processor attempts their trial in parallel,
  // without analyze modify or checkpoint.
  // Store new macrostate and acceptance prob
//  #pragma omp parallel private(proc_id)
//  {
//    proc_id = omp_get_thread_num();
  // serialize for debugging
  for (int proc_id = 0; proc_id < num_threads_; ++proc_id) {
    Pool * pool = &trial_pool_[proc_id];
    pool->set_accepted(clones_[proc_id].attempt_trial(pool->index(), random));
    pool->set_ln_prob(
      clones_[proc_id].trial(pool->index())->accept().ln_metropolis_prob()
    );
    INFO("p" << omp_get_thread_num() << " critical proc_id " << proc_id << " " << pool->str());
  }

  // Determine first that any trial was accepted.
  int first_thread_accepted = num_threads_;
  for (int ithread = 0; ithread < num_threads_; ++ithread) {
    if (trial_pool_[ithread].accepted() &&
        first_thread_accepted == num_threads_) {
      first_thread_accepted = ithread;
    }
  }
  INFO("first thread " << first_thread_accepted);

  // revert trials after accepted trial.
  if (first_thread_accepted != num_threads_) {
    ERROR("not tested");
    for (int ithread = first_thread_accepted + 1;
         ithread < num_threads_;
         ++ithread) {
      INFO("reverting trial " << ithread);
      const Pool& pool = trial_pool_[ithread];
      clones_[ithread].revert(pool.index(), pool.accepted());
    }
  }

  // for each thread up to the first accepted,
  // update other threads (incl. main) regarding failed attempt by thread.
  for (int ithread = 0; ithread < first_thread_accepted; ++ithread) {
    INFO("mimic failed " << ithread);
    const Pool& pool = trial_pool_[ithread];
    // loop through other threads to update
    for (int jthread = 0; jthread < num_threads_; ++jthread) {
      clone_(ithread, jthread)->mimic_trial_rejection(
        pool.index(), pool.ln_prob());
    }
  }

  // how to mimic a trial acceptance ?
  // record each random choice ?
  // use this record to perform the trial ?
  // set random seeds to be the same ?

  // finally, update all other threads regarding first accepted attempt
  if (first_thread_accepted < num_threads_) {
    ERROR("not implemented");
    for (int jthread = 0; jthread < num_threads_; ++jthread) {
//      clone_(first_thread_accepted, jthread)->mimic_trial_accepted(
//        clone_[first_thread_accepted]
//      );
      // HWH unoptimized but should work for testing
      *clone_(first_thread_accepted, jthread) = clones_[first_thread_accepted];
    }
  }

  // update analyze/modify/checkpoint for master after each update.
  after_trial_();
  for (int ithread = 0; ithread < num_threads_; ++ithread) {
    clones_[ithread].after_trial_();
  }

  // check that the current energy of all threads and master are the same.
  const double energy = criteria->current_energy();
  const double tolerance = 1e-8;
  const int num_particles = system->configuration().num_particles();
  for (int ithread = 0; ithread < num_threads_; ++ithread) {
    const double diff = clones_[ithread].criteria()->current_energy() - energy;
    ASSERT(fabs(diff) <= tolerance, "diff: " << diff);
    ASSERT(num_particles == clones_[ithread].system().configuration().num_particles(),
      "err");
  }

    // periodically equate all clones to main to prevent drift (energy checks?)
}



//    if (num_pool_ == -1) num_pool_ = num_threads_;
//    ASSERT(num_pool_ >= num_threads_,
//      "pool of trials(" << num_pool_ << ") must be greater than or equal to " <<
//      "number of threads " << num_threads_);
//  INFO(num_pool_);

//  INFO("num_pool_ " << num_pool_ << " num_threads_ " << num_threads_);
//  if (num_pool_ == num_threads_) {
//    distribute_for_();
//  } else {
//    distribute_while_();
//    ERROR("not implemented");
//  }

//void Pipeline::distribute_while_() {
//  int sj = -1;    // shared loop counter
//  int sstop = 0;  // shared loop condition
//  int tn, tj, tstop;
//
//  #pragma omp parallel private(tn, tj, tstop)
//  {
//    tn = omp_get_thread_num();
//    tstop = 0;
//    while (!sstop) {
//      #pragma omp critical
//      {
//        ++sj;     // increment shared loop counter
//        tj = sj;  // keep a private copy
//        if (tj == num_pool_ -1) {
//          sstop = 1;
//          #pragma omp flush(sstop)
//        }
//        INFO("critical tj " << tj << " tn " << tn);
//      }
//
//      // do work, update tstop to terminate
//      // clones_[tj].attempt_trial(trial_pool_[tj].index);
//
//      // check for termination
//      if (tstop) {
//        sstop = 1;
//        #pragma omp flush(sstop)
//      }
//    }
//  }
//  INFO("stopped");
//}

}  // namespace feasst
