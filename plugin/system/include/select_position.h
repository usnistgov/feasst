
#ifndef FEASST_SYSTEM_SELECT_POSITION_H_
#define FEASST_SYSTEM_SELECT_POSITION_H_

#include <vector>
#include "configuration/include/select.h"

namespace feasst {

// HWH move to configuration plugin?
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

  void add_site(const int particle_index, const int site_index) override {
    Select::add_site(particle_index, site_index);
    resize();
  }
//  /// Similar to Select::exchange_indices, but also positions and properties.
//  bool exchange_indices_positions(const SelectPosition& select) {
//    if (!exchange_indices(select)) return false;
//    for (int ipart = 0; ipart < static_cast<int>(site_positions_.size()); ++ipart) {
//      particle_positions_[ipart] = select.particle_positions()[ipart];
//      std::vector<Position> * sites = &site_positions_[ipart];
//      const std::vector<Position>& sel_sites = select.site_positions_[ipart];
//      if (sites->size() != sel_sites.size()) {
//        return false;
//      } else {
//        for (int isite = 0; isite < static_cast<int>(sites->size()); ++isite) {
//          // HWH opt consider implementing position/properties exchange
//          (*sites)[isite] = sel_sites[isite];
//          site_properties_[ipart][isite] = select.site_properties_[ipart][isite];
//        }
//      }
//    }
//    return true;
//  }

  // optimized access to site positions.
  Position * get_site_position(const int particle_index, const int site_index) {
    return &site_positions_[particle_index][site_index]; }

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
