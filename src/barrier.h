/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef BARRIER_H_
#define BARRIER_H_

#include "./base.h"
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

// Virtual class for different features of the barrier.
class Feature : public Base {
public:
	Feature () {}
	virtual ~Feature() {;}
};

// Overall, composite barrier composed of different features for a single particle type.
class Barrier : public Base {
public:
 	Barrier () {}
 	virtual ~Barrier() {}
};

// All barriers for each particle type.
class SpeciesBarriers : public Base {
public:
 	SpeciesBarriers () {}
 	~SpeciesBarriers() {;}
};

// "Hard" slit pore with no interactions with the species
class HardSlitPore : public Feature {
public:
 	HardSlitPore () {}
 	~HardSlitPore () {;}
};

// Slit pore with square well interactions
class SqwSlitPore : public Feature {
public:
	SqwSlitPore () {}
 	~SqwSlitPore () {;}
};

// Slit pore with square well interactions
class SqwCylinder : public Feature {
public:
	SqwCylinder () {}
 	~SqwCylinder () {;}
};

// Add additional barriers here ...

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // BARRIER_H_
