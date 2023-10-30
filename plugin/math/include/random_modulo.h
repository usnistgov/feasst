
#ifndef FEASST_MATH_RANDOM_MODULO_H_
#define FEASST_MATH_RANDOM_MODULO_H_

#include <memory>
#include <random>
#include <sstream>
#include "math/include/random.h"

namespace feasst {

/**
  This is a poor quality implementation of a random number generator,
  as described in Allen and Tildesley "Computer Simulation of Liquids"
  Appendix G.
  Because RandomMT19937 has known differences between GCC and Clang, this
  implementation allows exact reproduction on different compilers and
  operating systems.
  The higher-quality RandomMT19937 is recommended for production simulations.
 */
class RandomModulo : public Random {
 public:
  //@{
  /** @name Arguments
    - Random arguments.
   */
  explicit RandomModulo(argtype args = argtype());
  explicit RandomModulo(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  // serialize
  std::shared_ptr<Random> create(std::istream& istr) const override {
    return std::make_shared<RandomModulo>(istr); }
  std::shared_ptr<Random> create(argtype * args) const override {
    return std::make_shared<RandomModulo>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RandomModulo(std::istream& istr);
  virtual ~RandomModulo() {}
  //@}
 private:
  long int M_, a_, X_;

  void reseed_(const int seed) override;
  double gen_uniform_() override;
};

inline std::shared_ptr<RandomModulo> MakeRandomModulo(
    argtype args = argtype()) {
  return std::make_shared<RandomModulo>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_RANDOM_MODULO_H_
