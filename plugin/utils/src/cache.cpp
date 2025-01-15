#include "utils/include/cache.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"

namespace feasst {

Cache::Cache() {
  set_load();
  set_unload();
}

void Cache::set_load(const bool store) {
  is_loading_ = store;
  is_unloading_ = false;
  stored_.clear();
}

void Cache::set_unload(const bool unload) {
  is_unloading_ = unload;
  is_loading_ = false;
}

void Cache::set_unload(const Cache& cache) {
  set_unload(true);
  ASSERT(cache.is_loading_, "other cache was not storing values");
  //ASSERT(cache.stored_.size() > 0, "other cache has no stored values");
  stored_ = cache.stored_;
}

bool Cache::is_unloading(double * value) {
  if (is_unloading_) {
    ASSERT(stored_.size() > 0, "can not unload if nothing stored");
    *value = stored_.front();
    stored_.pop_front();
    return true;
  }
  return false;
}

void Cache::load(const double value) {
  if (is_loading_) {
    DEBUG("storing: " << value);
    stored_.push_back(value);
  }
}

void Cache::serialize(std::ostream& ostr) const {
  feasst_serialize_version(989, ostr);
  feasst_serialize(is_loading_, ostr);
  feasst_serialize(is_unloading_, ostr);
  feasst_serialize(stored_, ostr);
}

Cache::Cache(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 989, "mismatch version: " << version);
  feasst_deserialize(&is_loading_, istr);
  feasst_deserialize(&is_unloading_, istr);
  feasst_deserialize(&stored_, istr);
}

}  // namespace feasst
