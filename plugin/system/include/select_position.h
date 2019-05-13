
#ifndef FEASST_SYSTEM_SELECT_POSITION_H_
#define FEASST_SYSTEM_SELECT_POSITION_H_

#include <vector>
#include "configuration/include/select.h"

namespace feasst {

/**
  A selection which includes site and particle positions.
  Also include properties (e.g., for ewald, or neighbors, etc)
 */
class SelectPosition : public Select {
 public:
  SelectPosition() : Select() {}

  SelectPosition(const Select& select, const ParticleFactory& particles)
    : Select(select) {
    resize();
    load_positions(particles);
  }

  /// Return the site positions.
  const std::vector<std::vector<Position> >& site_positions() const {
    return site_positions_;
  }

  /// Return the site positions.
  const std::vector<std::vector<Properties> >& site_properties() const {
    return site_properties_;
  }

  /// Return the particle positions.
  const std::vector<Position>& particle_positions() const {
    return particle_positions_;
  }

  /// Set the position of a site by particle and site index.
  /// Note that these indices are based on selection, not configuration.
  void set_site_position(const int particle_index,
                         const int site_index,
                         const Position& position);

  /// Same as above except vector position is accepted.
  void set_site_position(const int particle_index,
                         const int site_index,
                         const std::vector<double> coord);

  /// Set the property of a site by particle and site index.
  /// Note that these indices are based on selection, not configuration.
  void set_site_properties(const int particle_index,
                           const int site_index,
                           const Properties& properties);

  /// Set the position of a particle by its index.
  /// Note that this index is based on selection, not configuration.
  void set_particle_position(const int particle_index,
                             const Position& position);

  /// Load the positions from the existing selection indices.
  void load_positions(const ParticleFactory& particles);

  void clear() override {
    Select::clear();
    clear_();
  }

  void resize();//const Select& select);

  /// Remove the last site.
  void remove_last_site() override;

  /// Remove the first site.
  void remove_first_site() override;

  virtual ~SelectPosition() {}

 private:
  std::vector<Position> particle_positions_;
  std::vector<std::vector<Position> > site_positions_;
  std::vector<std::vector<Properties> > site_properties_;

  void clear_() {
    particle_positions_.clear(); site_positions_.clear();
  }
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_SELECT_POSITION_H_
