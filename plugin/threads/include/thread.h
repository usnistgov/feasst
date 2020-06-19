#ifndef FEASST_THREAD_THREAD_H_
#define FEASST_THREAD_THREAD_H_

namespace feasst {

class Thread {
 public:
  Thread() {}

  /// Return true if parallel library is enabled.
  virtual bool is_enabled() const = 0;

  /// Return the number of threads.
  virtual int num() const = 0;

  /// Return the thread index;
  virtual int thread() const = 0;

  /// Return true if the current iteration of a large task that is broken into
  /// chunks is to be performed by this processor.
  bool in_chunk(const int iteration, const int total_iterations) const;

  virtual ~Thread() {}
};

}  // namespace feasst

#endif  // FEASST_THREAD_THREAD_H_
