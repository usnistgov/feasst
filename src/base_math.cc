/**
 * \file
 *
 * \brief
 */

#include "./base_math.h"

/**
 * Constructor
 */
BaseMath::BaseMath() {
  verbose_ = 0;
  className_.assign("BaseMath");
}

/**
 * random quaterion
 *  FRANZ J. VESELY, JOURNAL OF COMPUTATIONAL PHYSICS 47, 291-296 (1982)
 */
vector<double> BaseMath::quatRandom() {
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

/**
 * function to return random position within spherical shell defined by rabove
 * and rebelow in dim dimensions
 */
vector<double> BaseMath::ranShell(
  const double rabove,     //!< upper limit of shell
  const double rbelow,     //!< lower limit of shell
  const int dim            //!< spatial dimensions
  ) {
  vector<double> x(dim);

  // error check
  ASSERT(rabove >= 0, "rabove(" << rabove << ") must be positive");
  ASSERT(rbelow >= 0, "rbelow(" << rbelow << ") must be positive");
  ASSERT(rabove < 1000*rbelow, "not efficiently implemented for rabove ("
         << rabove << ") >> (" << rbelow << ")");
  ASSERT(rabove >= rbelow, "rabove (" << rabove << ") must be greater than "
         << "rebelow (" << rbelow << ")");
  x = ranUnitSphere(3);
  const double fac = pow(uniformRanNum() * (rabove*rabove*rabove
    - rbelow*rbelow*rbelow) + rbelow*rbelow*rbelow, 1./3.);
  for (vector<double>::iterator it = x.begin(); it != x.end(); ++it) {
    *it *= fac;
  }
  return x;
}

/**
 * given cumulative discrete probability distribution in vector of doubles,
 * return chosen integer element from uniform probability distribution
 */
int BaseMath::ranFromCPDF(
  const vector<double> &cpdf  //!< cumulative probability distriubtion function
  ) {
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

/**
 * generate a randomly preturbed quaterion
 */
vector<double> BaseMath::quatRandom(
  const double maxDisp) { //!< amount to preturb quaterion from unit vector
                          // (e.g. identity rotation matrix)
  vector<double> q = quatRandom();
  vector<double> qunit(4);
  qunit[3] = 1.;
  for (unsigned int i = 0; i < q.size(); ++i) {
    q[i] = qunit[i] + maxDisp*q[i];
  }
  const double qsize = sqrt(myVecDotProd(q, q));
  for (unsigned int i = 0; i < q.size(); ++i) {
    q[i] /= qsize;
  }
  return q;
}

/**
 * function to return random position on unit sphere in dim dimensions
 */
vector<double> BaseMath::ranUnitSphere(
  const int dim) {  //!< spatial dimensions
  ASSERT(dim == 3,
    "only implemented for 3 dimensions (dim=" << dim << ") using quaterions");
  vector<double> q = quatRandom();
  vector<vector<double> > r = quat2rot(q);
  vector<double> a;
  a.push_back(1); a.push_back(0); a.push_back(0);
  vector<double> x = myMatVecMul(r, a);
  return x;
}

/**
 * return angle selected from probability distribution associated with
 * harmonic bond bending energy, \beta*U=k0*(t-t0)**2
 * Frenkel and Smit, page 343, below Equation 13.3.6
 */
double BaseMath::ranAngle(const double k0, const double t0, const int power) {
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

/**
 * return bond length selected from probability distribution associated with
 * harmonic bond bending energy, \beta*U=k0*(l-l0)
 */
double BaseMath::ranBond(const double k0, const double l0, const int power) {
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

/**
 * return random euler angles which yeild uniform sampling on a surface of a
 * sphere (e.g., uniform random orientation)
 *
 * Euler angles (phi, theta, psi) are defined by rotation phi about z-axis,
 * rotation theta about new x-xis, then rotation psi about new z-axis
 *  phi [-pi,pi], theta [0, pi], psi [-pi, pi]
 */
vector<double> BaseMath::eulerRandom() {
  vector<double> euler(3);
  euler[0] = PI*(2.*uniformRanNum() - 1.);
  euler[1] = acos(uniformRanNum());
  euler[2] = PI*(2.*uniformRanNum() - 1.);
  return euler;
}

