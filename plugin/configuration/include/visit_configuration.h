
#ifndef FEASST_SYSTEM_VISIT_CONFIGURATION_H_
#define FEASST_SYSTEM_VISIT_CONFIGURATION_H_

#include "configuration/include/visit_particles.h"

namespace feasst {

class Configuration;

struct LoopDescriptor {
  int particle_index;
  int site_index;
  int particle_type;
};

class LoopConfigOneBody;

class VisitConfiguration : public VisitParticles {
 public:
  void loop(const Configuration& config,
            LoopConfigOneBody * loop_config_one_body,
            const Select& select);
  void loop(const Configuration& config,
            LoopConfigOneBody * loop_config_one_body,
            const int group_index = 0);

  virtual ~VisitConfiguration() {}
 private:
  LoopDescriptor data_;
};

class LoopConfigOneBody {
 public:
  void visit(VisitConfiguration * visitor,
             const Configuration& configuration,
             const Select& select) {
    visitor->loop(configuration, this, select);
  }
  virtual void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) = 0;
  virtual ~LoopConfigOneBody() {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_CONFIGURATION_H_
