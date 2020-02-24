
#ifndef FEASST_CHAIN_UTILS_CHAIN_H_
#define FEASST_CHAIN_UTILS_CHAIN_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
#include "chain/include/trial_deprotonation.h"
#include "chain/include/trial_protonation.h"

namespace feasst {

/**
  Add both deprotonation and protonation reactions with the same weight
  to ensure detailed balance.
  Arguments should match those documented for deprotonation.
  Protonation reaction arguments are then generated in this function.
 */
inline void add_deprotonation(MonteCarlo * monte_carlo,
    const argtype& args = argtype()) {
  monte_carlo->add(MakeTrialDeprotonation(args));

  // add protonation trial after generating the new args
  argtype prot_args = args;
  // first, rename add_type to remove_type
  std::map<std::string, std::string>::iterator it =
    prot_args.find("add_type");
//  it->second = "remove_type";
  prot_args.insert(std::pair<std::string, std::string>("remove_type",
                                                         it->second));
  prot_args.erase(it);

  // second, swap the reactant_site_type and new_site_type arguments
  Arguments args_(args);
  args_.dont_check();
  const std::string rst = args_.key("reactant_site_type").str();
  const std::string nst = args_.key("new_site_type").str();
  it = prot_args.find("reactant_site_type");
  prot_args.erase(it);
  it = prot_args.find("new_site_type");
  prot_args.erase(it);
  prot_args.insert(std::pair<std::string, std::string>("new_site_type", rst));
  prot_args.insert(std::pair<std::string, std::string>("reactant_site_type", nst));
  monte_carlo->add(MakeTrialProtonation(prot_args));
}

}  // namespace feasst

#endif  // FEASST_CHAIN_UTILS_CHAIN_H_
