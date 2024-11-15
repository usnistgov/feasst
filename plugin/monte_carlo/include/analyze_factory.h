
#ifndef FEASST_MONTE_CARLO_ANALYZE_FACTORY_H_
#define FEASST_MONTE_CARLO_ANALYZE_FACTORY_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

class TimerRDTSC;

/**
  Contains multiple Analyze objects.
 */
class AnalyzeFactory : public Analyze {
 public:
  explicit AnalyzeFactory(argtype args = argtype());

  void initialize(MonteCarlo * mc) override;

  /// Add an Analyze object.
  void add(std::shared_ptr<Analyze> analyze);

  /// Remove an Analyze object.
  void remove(const int index);

  /// Return the number.
  int num() const { return static_cast<int>(analyzers_.size()); }

  /// Return the Analyze objects.
  const std::vector<std::shared_ptr<Analyze> >& analyzers() const override {
    return analyzers_; }

  /// Return an Analyze object by index.
  const Analyze& analyze(const int index) const override {
    return const_cast<Analyze&>(*analyzers_[index]); }

  /// Write all Analyze immediately.
  void write_to_file(const MonteCarlo& mc) override;

  /// For use with CollectionMatrixSplice, transfer multistate between threads.
  void adjust_bounds(const bool adjusted_up, const std::vector<int>& states,
    AnalyzeFactory * analyze_factory);

  void trial(const MonteCarlo& mc) override;

  Analyze * get_analyze(const int index) override {
    return analyzers_[index].get(); }

  /// Set the timer
  void set_timer();

  /// Return timer
  const TimerRDTSC * const timer() const { return timer_.get(); }

  // serialize
  std::string class_name() const override {
    return std::string("AnalyzeFactory"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeFactory>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit AnalyzeFactory(std::istream& istr);
  virtual ~AnalyzeFactory();

 private:
  std::vector<std::shared_ptr<Analyze> > analyzers_;
  std::unique_ptr<TimerRDTSC> timer_;

  void trial_(const MonteCarlo& mc, const int index);

  int min_block_(const Criteria& criteria) const;
  std::string write_blocks_(const int min_block, const Accumulator& acc) const;
};

inline std::shared_ptr<AnalyzeFactory> MakeAnalyzeFactory(
    argtype args = argtype()) {
  return std::make_shared<AnalyzeFactory>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ANALYZE_FACTORY_H_
