import feasst

def system(config=None, box_length=8.109613, alphaL = 6.870983963962610000, kmax_squared=38, rcut=4.891304347826090):
    if not config:
      config = feasst.Configuration(feasst.MakeDomain(feasst.args({"cubic_box_length": str(box_length)})),
        feasst.args({"particle_type0": feasst.install_dir() + "/plugin/ewald/forcefield/data.rpm_plus",
                     "particle_type1": feasst.install_dir() + "/plugin/ewald/forcefield/data.rpm_minus"}))
    config.set_model_param("cutoff", 0, rcut)
    config.set_model_param("cutoff", 1, rcut)
    system = feasst.System()
    system.add(config)
    system.add(feasst.Potential(feasst.MakeEwald(feasst.args({"kmax_squared": str(kmax_squared),
        "alpha": str(alphaL/system.configuration().domain().min_side_length())}))))
# Unfortunatley, swig isn't accepting the below method of constructing a two body factory
#    system.add(feasst.Potential(feasst.MakeModelTwoBodyFactory(
#        feasst.ModelTwoBodyVector([feasst.MakeLennardJones(), feasst.MakeChargeScreened()]))))
    two = feasst.MakeModelTwoBodyFactory()
    two.add(feasst.MakeHardSphere())
    two.add(feasst.MakeChargeScreened())
    system.add(feasst.Potential(two))
    system.add(feasst.Potential(feasst.MakeChargeSelf()))
    # system.precompute()
    return system
