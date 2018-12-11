
/**

Bonds # owned by particle (type + sites)
ModelBonds # owned by system
#BondParams # owned by particles or config (like model params)
## or not bondparams.. put in custom properties of bonds.. pass bonds

bond parameters -> spring constant, equilibrium distance/angle, etc

visit_bonds(config, ModelBondFactory):

#optimize: find indices of parameters needed for each bond model/property combo

for each particle
  type = particle.type
  for each bond in particle type
    bondModel = ModelBondFactory(bondtype)
    bond = particle(type).bond
    sites = particle.sites(bond) # list of constant pointers/reference?
    bondParams = bond.params (or model params)
    bondModel->compute(bondParameters, sites)

    bondModel::compute(...)
      ASSERT(sites.size() == ...)
      lx = (x-d...)
 */
