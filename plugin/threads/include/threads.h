#ifndef FEASST_THREADS_THREADS_H_
#define FEASST_THREADS_THREADS_H_

namespace feasst {

class Threads {
 public:
  /// Return true if parallel library is enabled.
  virtual bool is_enabled() const = 0;

  /// Return the number of threads.
  virtual int num() const = 0;

  virtual ~Threads() {}
};

}  // namespace feasst

#endif  // FEASST_THREADS_THREADS_H_
