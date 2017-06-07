#include "./base_random.h"
#include "./random_nr3.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

BaseRandom::BaseRandom() {
  verbose_ = 0;
  className_.assign("BaseRandom");
  clearRNG();
}

BaseRandom::BaseRandom(const char* fileName) {
  verbose_ = 0;
  className_.assign("BaseRandom");
  initRNG(fileName);
}

void BaseRandom::writeRngRestart(const char* fileName) {
  if (ranNum_) {
    stringstream ss;
    ss << fileName << "rng";
    ranNum_->writeRestart(ss.str().c_str());
  }
}

void BaseRandom::initRNG(shared_ptr<Random> ran) {
  ranNum_ = ran;
}

void BaseRandom::initRNG(unsigned long long seed) {
  if (seed == 0) seed = rand();
  ranNum_ = std::make_shared<RandomNR3>(seed);
}

void BaseRandom::initRNG(const char* fileName) {
  stringstream ss;
  ss << fileName << "rng";
  if (fileExists(ss.str().c_str())) {
    ranNum_ = std::make_shared<RandomNR3>(ss.str().c_str());
  } else {
    clearRNG();
  }
}

void BaseRandom::clearRNG() {
  std::shared_ptr<Random> emptyRanNum;
  initRNG(emptyRanNum);
}

double BaseRandom::uniformRanNum() {
  if (!ranNum_) initRNG();
  const double ran = ranNum_->uniform();
  WARN(ran == 0., "uniformRanNum is exactly 0 to double precision accuracy");
  return ran;
}

int BaseRandom::uniformRanNum(const int min, const int max) {
  if (!ranNum_) initRNG();
  return ranNum_->uniform(min, max);
}

double BaseRandom::stdNormRanNum() {
  const double u = uniformRanNum();
  const double v = uniformRanNum();
  return sqrt(-2.*log(u))*cos(2.*PI*v);
}

vector<double> BaseRandom::stdRanNum(const double variance, const int dimen) {
  vector<double> ran(dimen);
  for (int dim = 0; dim < dimen; ++dim) {
    ran[dim] = variance*stdNormRanNum();
  }
  return ran;
}

double BaseRandom::gaussRanNumFS() {
  double r = 2., v1 = 0, v2 = 0;
  while (r >= 1) {
    v1 = 2.*uniformRanNum() - 1;
    v2 = 2.*uniformRanNum() - 1;
    r = v1*v1+v2*v2;
  }
  return v1*sqrt(-2.*log(r)/r);
}

// http://stackoverflow.com/questions/440133/
// how-do-i-create-a-random-alpha-numeric-string-in-c
std::string BaseRandom::randomHash(const int length) {
  std::string str;
  const std::string alphanum = "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for (int i = 0; i < length; ++i) {
    std::stringstream ss;
    ss << alphanum[uniformRanNum(0, alphanum.size()-1)];
    str.append(ss.str());
  }
  return str;
}

vector<double> BaseRandom::quatRandom() {
  vector<double> q(4);
  double s1 = 2;
  double s2 = 2;
  while (s1 > 1) {
    q[0] = uniformRanNum();
    q[0] = 2*q[0] - 1;
    q[1] = uniformRanNum();
    q[1] = 2*q[1] - 1;
    s1 = q[0]*q[0] + q[1]*q[1];
  }
  while (s2 > 1) {
    q[2] = uniformRanNum();
    q[2] = 2*q[2] - 1;
    q[3] = uniformRanNum();
    q[3] = 2*q[3] - 1;
    s2 = q[2]*q[2] + q[3]*q[3];
  }
  s1 = sqrt((1 - s1)/s2);
  q[2] = q[2]*s1;
  q[3] = q[3]*s1;
  return q;
}

vector<double> BaseRandom::quatRandom(
  const double maxPerturb) {
  vector<double> q = quatRandom();
  vector<double> qunit(4);
  qunit[3] = 1.;
  for (unsigned int i = 0; i < q.size(); ++i) {
    q[i] = qunit[i] + maxPerturb*q[i];
  }
  const double qsize = sqrt(vecDotProd(q, q));
  for (unsigned int i = 0; i < q.size(); ++i) {
    q[i] /= qsize;
  }
  return q;
}

vector<double> BaseRandom::ranShell(const double rabove, const double rbelow,
  const int dim) {
  vector<double> x(dim);
  ASSERT(rabove >= 0, "rabove(" << rabove << ") must be positive");
  ASSERT(rbelow >= 0, "rbelow(" << rbelow << ") must be positive");
  ASSERT(rabove < 1000*rbelow, "not efficiently implemented for rabove ("
         << rabove << ") >> (" << rbelow << ")");
  ASSERT(rabove >= rbelow, "rabove (" << rabove << ") must be greater than "
         << "rebelow (" << rbelow << ")");
  ASSERT(dim == 3, "ranShell implemented for 3D only.");
  x = ranUnitSphere(3);
  const double fac = pow(uniformRanNum() * (rabove*rabove*rabove
    - rbelow*rbelow*rbelow) + rbelow*rbelow*rbelow, 1./3.);
  for (vector<double>::iterator it = x.begin(); it != x.end(); ++it) {
    *it *= fac;
  }
  return x;
}

vector<double> BaseRandom::ranUnitSphere(const int dim) {
  ASSERT(dim == 3,
    "only implemented for 3 dimensions (dim=" << dim << ") using quaterions");
  vector<double> q = quatRandom();
  vector<vector<double> > r = quat2rot(q);
  vector<double> a;
  a.push_back(1); a.push_back(0); a.push_back(0);
  vector<double> x = matVecMul(r, a);
  return x;
}

int BaseRandom::ranFromCPDF(const vector<double> &cpdf) {
  double prev = -1, current = -1;
  const double ranNum = uniformRanNum();
  int selected = 0, i = -1;
  ASSERT(fabs(cpdf.back()-1) < doubleTolerance, "cpdf must end in 1 ("
         << cpdf.back() << " != 1)");
  while (selected == 0) {
    ++i;
    current = cpdf[i];
    if (ranNum < current) {
      selected = 1;
    }
    ASSERT(current >= prev, "cpdf must be monotonically nondecreasing");
    ASSERT((i+1 != static_cast<int>(cpdf.size())) || (selected != 0),
      "no selection in cpdf"
      << std::endl << i << " " << current << " " << prev << " " << cpdf.back());
    prev = current;
  }
  return i;
}

double BaseRandom::ranAngle(const double k0, const double t0, const int power) {
  int accept = 0, i = 0;
  double theta = 0;
  const int maxAttempts = 1e6;
  double minAngle = 0.;
  if (fabs(k0) < 100*doubleTolerance) minAngle = t0;
  while ( (accept == 0) && (i < maxAttempts) ) {
    theta = minAngle + (PI-minAngle)*uniformRanNum();
    const double dtheta = theta-t0;
    if (uniformRanNum() < sin(theta)*exp(-k0*pow(dtheta, power))) {
      accept = 1;
    }
    ++i;
  }
  ASSERT(i != maxAttempts,
    "maximum attempts reached in selected ranAngle(" << i << ")");
  return theta;
}

double BaseRandom::ranBond(const double k0, const double l0, const int power) {
  int accept = 0, i = 0;
  double length = 0;
  const int maxAttempts = 1e6;
  double lmax = l0*2.;
  const double lmin = 0.;
  if (fabs(k0) < 100*doubleTolerance) lmax = l0;

  // alternative 1, uniform and then accept/reject
  while ( (accept == 0) && (i < maxAttempts) ) {
    length = lmin + (lmax-lmin)*uniformRanNum();
    const double dl = length-l0;
    if (uniformRanNum() < pow(length/lmax, 2)*exp(-0.5*k0*pow(dl, power))) {
      accept = 1;
    }

//  // alternative 2, direct guassian if power == 2
//  const double sigma = sqrt(1./k0);
//  const double a = pow(l0+3.*sigma, 2);
//  ASSERT(power == 2, "assumes gaussian");
//  while ( (accept == 0) && (i < maxAttempts) ) {
//    length = gaussRanNum(sigma, l0);
//    if ( (length >= lmin) && (length <= lmax) ) {
//      if (uniformRanNum() < length*length/a) accept = 1;
//    }

  // end alternatives
    ++i;
  }
  ASSERT(i != maxAttempts,
    "maximum attempts reached in selected ranBond(" << i << ")");
  return length;
}

vector<double> BaseRandom::eulerRandom() {
  vector<double> euler(3);
  euler[0] = PI*(2.*uniformRanNum() - 1.);
  euler[1] = acos(uniformRanNum());
  euler[2] = PI*(2.*uniformRanNum() - 1.);
  return euler;
}

void BaseRandom::reconstructDerived_() {
  clearRNG();
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


