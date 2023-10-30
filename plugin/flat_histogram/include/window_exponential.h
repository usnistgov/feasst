
#ifndef FEASST_FLAT_HISTOGRAM_WINDOW_EXPONENTIAL_H_
#define FEASST_FLAT_HISTOGRAM_WINDOW_EXPONENTIAL_H_

#include <vector>
#include "utils/include/arguments.h"
#include "flat_histogram/include/window.h"

namespace feasst {

/**
  Determine windows based on an expontential factor, \f$\alpha\f$.
  The macrostate of the windows, \f$m_i\f$, is given by

\f$m_{i+1} = \left( m_i^{\alpha} + \frac{m_{n+1}^{\alpha} - m_0^{\alpha}}{n} \right)^{1/\alpha}\f$

  where \f$i\f$ is the index of the window, ranging from \f$0\f$ to \f$n+1\f$
  and \f$n\f$ is the number of windows.
  Thus, \f$m_{n+1}\f$ is the maximum macrostate over the entire range, while
  \f$m_0\f$ is the minimum macrostate over the entire range.
  For \f$\alpha=1\f$, the widths of all windows are the same.
  For \f$\alpha>1\f$, the widths of the windows decrease exponentially.
  For \f$\alpha<1\f$, the widths of the windows increase exponentially.

  The given formula is continuous and not rounded.
  Instead, use Window::boundaries for handling rounding and overlap.
  For example, 4 windows over a macrostate range of [0, 200], \f$\alpha=2\f$,
  results in the following segments:
  [0, 100.00, 141.42, 173.21, 200].

  The choice of alpha can affect the efficiency of the simulation, depending on
  the relative speed of the convergence of the low macrostates compared to the
  higher macrostates.
  If the lower macrostate windows are converging slower, then alpha should be
  decreased.
  If the higher macrostate simulations are converging slower, then alpha may
  need to be increased.
  The choice of alpha is dependent on the system, conditions and convergence
  criteria.
 */
class WindowExponential : public Window {
 public:
  //@{
  /** @name Arguments
    - alpha: exponential factor (default: 1.5).
    - min[i]: minimum macrostate value of the i-th window (default: None),
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      This allows the user to override alpha at the smaller windows.
    - Window arguments.
   */
  explicit WindowExponential(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  std::vector<double> segment() const override;
  virtual ~WindowExponential() {}

  //@}
 private:
  double alpha_;
  std::vector<double> partial_segment_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WINDOW_EXPONENTIAL_H_
