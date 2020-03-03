import feasst

def system(config, alphaL = 5.6, kmax_squared=38, rcut=False):
    config.add_model_param("alpha", alphaL/config.domain().min_side_length())
    if rcut:
        config.set_model_param("cutoff", 0, rcut)
        config.set_model_param("cutoff", 1, rcut)
    system = feasst.System()
    system.add(config)
    system.add(feasst.Potential(feasst.MakeEwald(feasst.args({"kmax_squared": str(kmax_squared),
        "alpha": str(5.6/system.configuration().domain().min_side_length())}))))
# Unfortunatley, swig isn't accepting the below method of constructing a two body factory
#    system.add(feasst.Potential(feasst.MakeModelTwoBodyFactory(
#        feasst.ModelTwoBodyVector([feasst.MakeLennardJones(), feasst.MakeChargeScreened()]))))
    two = feasst.MakeModelTwoBodyFactory()
    two.add(feasst.MakeLennardJones())
    two.add(feasst.MakeChargeScreened())
    system.add(feasst.Potential(two))
    system.add(feasst.Potential(feasst.MakeChargeScreenedIntra(), feasst.MakeVisitModelBond()));
    system.add(feasst.Potential(feasst.MakeChargeSelf()));
    system.add(feasst.Potential(feasst.MakeLongRangeCorrections()));
    # system.precompute()
    return system
