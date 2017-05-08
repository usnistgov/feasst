/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 */

#include "./base_random.h"
#include "./random_nr3.h"

namespace feasst {

/**
 * Constructor
 */
BaseRandom::BaseRandom() {
  verbose_ = 0;
  className_.assign("BaseRandom");
  ranNumOwned_ = false;
  clearRNG();
}
BaseRandom::BaseRandom(const char* fileName) {
  verbose_ = 0;
  className_.assign("BaseRandom");
  ranNumOwned_ = false;
  initRNG(fileName);
}

BaseRandom::~BaseRandom() {
  if (ranNumOwned_) delete ranNum_;
}

/**
 * reset object pointers
 */
void BaseRandom::reconstructDerived() {
  ranNumOwned_ = false;
  clearRNG();
}

/**
 * initialize random number generator
 */
void BaseRandom::initRNG(Random *ran) {
  if (ranNumOwned_) delete ranNum_;
  ranNum_ = ran;
  ranNumOwned_ = false;
}

/**
 * initialize random number generator
 *  seed with compiler random number generator
 */
void BaseRandom::initRNG(unsigned long long seed) {
  if (ranNumOwned_) delete ranNum_;
  if (seed == 0) seed = rand();
  ranNum_ = new RandomNR3(seed);
  ranNumOwned_ = true;
}

/**
 * initialize random number generator
 */
void BaseRandom::initRNG() {
  // zero seed signals use of compiler random number for seed
  initRNG(static_cast<unsigned long long>(0));
}

/**
 * clear random number generator
 *  typically used after clone because no deletion of pointer memory
 */
void BaseRandom::clearRNG() {
  initRNG(reinterpret_cast<Random*>(NULL));
}

/**
 * initialize random number generator from file
 */
void BaseRandom::initRNG(const char* fileName) {
  stringstream ss;
  ss << fileName << "rng";
  if (fileExists(ss.str().c_str())) {
    if (ranNumOwned_) delete ranNum_;
    ranNum_ = new RandomNR3(ss.str().c_str());
    ranNumOwned_ = true;
  } else {
    clearRNG();
  }
}

/**
 * uniform random doubleprecision number between 0 and 1
 */
double BaseRandom::uniformRanNum() {
  if (ranNum_ == NULL) initRNG();
  const double ran = ranNum_->uniform();
  WARN(ran == 0., "uniformRanNum is exactly 0 to double precision accuracy");
  return ran;
}

/**
 * uniform random integer between range min and max, inclusive
 */
int BaseRandom::uniformRanNum(const int min, const int max) {
  if (ranNum_ == NULL) initRNG();
  return ranNum_->uniform(min, max);
}

/**
 * standard normal random doubleprecision number
 *  uses Box-Muller method
 */
double BaseRandom::stdNormRanNum() {
  const double u = uniformRanNum();
  const double v = uniformRanNum();
  return sqrt(-2.*log(u))*cos(2.*PI*v);
}

/**
 * gaussian random number of given variance
 */
vector<double> BaseRandom::stdRanNum(const double variance,
  const int dimen     //!< dimension
  ) {
  vector<double> ran(dimen);
  for (int dim = 0; dim < dimen; ++dim) {
    ran[dim] = variance*stdNormRanNum();
  }
  return ran;
}

/**
 * write restart file
 */
void BaseRandom::writeRngRestart(const char* fileName) {
  if (ranNum_ != NULL) {
    stringstream ss;
    ss << fileName << "rng";
    ranNum_->writeRestart(ss.str().c_str());
  }
}

/// gassuian distribution random number from Frenkel and Smit
double BaseRandom::gaussRanNumFS() {
  double r = 2., v1 = 0, v2 = 0;
  while (r >= 1) {
    v1 = 2.*uniformRanNum() - 1;
    v2 = 2.*uniformRanNum() - 1;
    r = v1*v1+v2*v2;
  }
  return v1*sqrt(-2.*log(r)/r);
}

/// generate a random hash
std::string BaseRandom::randomHash(const int len) {
  std::string str;
  const std::string alphanum = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for (int i = 0; i < len; ++i) {
    std::stringstream ss;
    ss << alphanum[uniformRanNum(0, alphanum.size()-1)];
    str.append(ss.str());
  }
  return str;
}

}  // namespace feasst

