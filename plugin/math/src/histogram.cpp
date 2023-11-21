
#include <cmath>
#include "math/include/histogram.h"
#include "utils/include/debug.h"
#include "math/include/formula_polynomial.h"
#include "math/include/utils_math.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"  // NEAR_ZERO
#include "math/include/formula.h"

namespace feasst {

Histogram::Histogram(argtype args) : Histogram(&args) { FEASST_CHECK_ALL_USED(args); }
Histogram::Histogram(argtype * args) {
  // optionally construct a constant width bin
  if (used("width", *args)) {
    const double width = dble("width", args);
    DEBUG("width " << MAX_PRECISION << width);
    const double max = dble("max", args);
    const double min = dble("min", args, 0.);
    ASSERT(max >= min, "max(" << max <<") should be >= min(" << min << ")");
    set_width_min_max(width, min, max);
  }
}

void Histogram::set_bin_size(const std::shared_ptr<Formula> bin_size) {
  set_expandable_();
  bin_size_ = bin_size;
  ASSERT(edges_.size() == 0 && size() == 0, "formula cannot be changed");
  histogram_.push_back(0.);
  const double fabove = bin_size_->evaluate(1);
  const double f0 = bin_size_->evaluate(0);
  const double fbelow = bin_size_->evaluate(-1);
  edges_.push_back(0.5*(fbelow + f0));
  edges_.push_back(0.5*(fabove + f0));
}

void Histogram::set_width_center(const double width, const double center) {
  auto bin_size = MakeFormulaPolynomial({{"x0", "0"}});
  bin_size->set_A(0, center).set_A(1, width);
  set_bin_size(bin_size);
  is_constant_width_ = 1;
}

void Histogram::set_width_min_max(const double width, const double min,
    const double max) {
  set_width_center(width, min);
  const double num = (max - min)/width;
  ASSERT(std::abs(num - round((num))) < 1e-8,
    "min: " << MAX_PRECISION << min << " and max: "
    << max << " do not align with width: " << width);
  for (double val = min; val < max + 0.5*width; val += width) {
    add(val, false);
  }
}

void Histogram::set_edges(const std::deque<double> edges) {
  ASSERT(edges.size() > 1, "must be more than 1 edge");
  edges_ = edges;
  set_not_expandable();
  histogram_.resize(edges_.size() - 1);
}

void Histogram::set_edges(const std::vector<double> edges) {
  std::deque<double> deque_edges;
  for (double edge : edges) deque_edges.push_back(edge);
  set_edges(deque_edges);
}

int Histogram::bin(const double value) const {
  if (is_constant_width_ == 1) {
    const double kWidth = edges_[1] - edges_[0];
    TRACE("kWidth " << kWidth);
    TRACE("value " << value);
    TRACE("center " << center_of_bin(0));
    return round((value - center_of_bin(0))/kWidth);
  } else {
    ASSERT(value <= max() && value >= min(),
      "histogram value(" << value << ") is out of range. max(" << max() <<
      ") min( " << min() << ")");
  }
  int bin = 0;
  bool found = false;
  while (!found && bin < size()) {
    if (value < edges_[bin + 1]) {
      found = true;
    }
    ++bin;
  }
  ASSERT(found == true, "out of range");
  return bin;
}

double Histogram::center_of_bin(const int bin) const {
  ASSERT(bin < size(), "bin: " << bin << " >= size: " << size());
  ASSERT(bin < static_cast<int>(edges_.size()),
    "bin: " << bin << " >= size: " << edges_.size());
  ASSERT(bin >= 0, "bin: " << bin << " < 0");
  return 0.5*(edges_[bin] + edges_[bin + 1]);
}

void Histogram::add(const double value, const bool update) {
  // initialize histogram if not already and formula is set
  ASSERT(edges_.size() != 0, "size error");
  if ((value <= max()) && (value >= min())) {
    if (update) ++histogram_[bin(value)];
  } else {
    // expand histogram if allowed and formula is available.
    ASSERT(bin_size_ != NULL && expandable_, "out of range");
    int increment = 0;
    const int kMaxIncrements = 1e8;
    int found = 0;
    if (value > max()) {
      while (found == 0 && increment < kMaxIncrements) {
        const double new_width = 2*(bin_size_->evaluate(size()) -
                                    edges_.back());
        edges_.push_back(edges_.back() + new_width);
        if (value < max()) {
          found = 1;
          if (update) {
            histogram_.push_back(1);
          } else {
            histogram_.push_back(0);
          }
        } else {
          histogram_.push_back(0);
        }
        ++increment;
      }
      ASSERT(increment != kMaxIncrements, "size error");
    } else if (value < min()) {
      while (found == 0 && increment < kMaxIncrements) {
        const double new_width = 2*(edges_.front() - bin_size_->evaluate(-1));
        edges_.push_front(edges_.front() - new_width);
        if (value > min()) {
          found = 1;
          if (update) {
            histogram_.push_front(1);
          } else {
            histogram_.push_front(0);
          }
        } else {
          histogram_.push_front(0);
        }
        bin_size_->set_x0(bin_size_->x0() + 1);
        ++increment;
      }
      ASSERT(increment != kMaxIncrements, "size error");
    } else {
      ASSERT(increment != -1, "strange range");
    }
  }
}

void Histogram::serialize(std::ostream& ostr) const {
  feasst_serialize_version(4226, ostr);
  feasst_serialize(histogram_, ostr);
  feasst_serialize(edges_, ostr);
  feasst_serialize(expandable_, ostr);
  feasst_serialize_fstdr(bin_size_, ostr);
  feasst_serialize(is_constant_width_, ostr);
}

Histogram::Histogram(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4226, "unrecognized verison: " << version);
  feasst_deserialize(&histogram_, istr);
  feasst_deserialize(&edges_, istr);
  feasst_deserialize(&expandable_, istr);
  // feasst_deserialize_fstdr(bin_size_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      bin_size_ = bin_size_->deserialize(istr);
    }
  }
  feasst_deserialize(&is_constant_width_, istr);
}

const std::string Histogram::str() const {
  std::stringstream ss;
  ss << "bin,num" << std::endl;
  for (int bin = 0; bin < size(); ++bin) {
    ss << center_of_bin(bin) << "," << histogram()[bin] << std::endl;
  }
  return ss.str();
}

}  // namespace feasst
