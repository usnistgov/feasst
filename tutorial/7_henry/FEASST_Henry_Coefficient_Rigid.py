import math
import feasst
import pyfeasst
import numpy as np
import multiprocessing as mp
import json

def process(in_file_name):
    with open(in_file_name + ".xyz", "r") as f:
        content = f.readlines()

    out_file_name = in_file_name + "_out"
    out_file = open(out_file_name + ".xyz", "w")
    print(content[0], file=out_file, end='')
    print("-1 " + content[1], file=out_file, end='')

    data = {}
    atomTypes = list()
    #numTypes = list()
    for line in range(len(content)):
        if line >= 2:
            name = content[line].lstrip().split(" ")[0]
            #print("name", name)
            if name not in atomTypes:
                atomTypes.append(name)
                data[name] = 1
            else:
                data[name] = data[name] + 1
    #print(data)

    for atom in atomTypes:
        for line in range(len(content)):
            if line >= 2:
                name = content[line].lstrip().split(" ")[0]
                if name == atom:
                    print(content[line], file=out_file, end='')

    with open(out_file_name + "_types.json", 'w') as json_file:
        json.dump(data, json_file)
    return data

def HenryCoefficient_worker(debug=False,seed=123456,**kwargs):

#     print("******************************* NOTE *******************************")
#     print("* Beware the subtle difference in xyz file formats.")
#     print("* FEASST uses the second line as [order_param, lx, ly, lz]")
#     print("* Dan uses the second line as [lx, ly, lz]")
#     print("* MFI_replicate.xyz follows Dans format,")
#     print("* while MFI_replicate_out.xyz from process.py follows FEASSTs format")
#     print("* co2.xyz and the input/outputs of concatenate.py follow FEASSTs format")
#     print("******************************* NOTE *******************************")

    if debug:
        input_parameters = { "adsorptive": "data.Ar",
                         "adsorbent": "./MFI_replicate",
                         "temperature": 350.,
                         "rcut": 15.,
                         "ncoeffs": 5,
                         "trials": 1.e4,
                         "pair_type": "LJCoulEwald",
                         "tail_type": "LFS",
                         "Ewald":{ "k2max": 27, "alpha": 5.6},
                         "scale": {"active": False, "factor": 10.},
                         "progress_bar": True
                           }
        feasst.ranInitForRepro( seed )
    else:
        # Reassmble the input kwargs into a dictionary
        input_parameters = dict(kwargs.items())
        try:
            feasst.ranInitForRepro( seed )
        except:
            feasst.ranInitByDate()

    #----------------------------
    #  MODEL PARAMETERS
    temp = input_parameters["temperature"]   # kelvin, assumes epsilons are kJ/mol
    adsorbate_filename = input_parameters["adsorptive"]
    activ = math.exp(-3.) # activity of adsorptive doesn't matter, but error checked
    rcut = input_parameters["rcut"]
    #-----------------------------

    space = feasst.makeSpace(feasst.args(
        {"dimen" : "3"}))
    in_file_name = input_parameters["adsorbent"]

    # potential Type
    if input_parameters["pair_type"] == "LJCoulEwald":
        pair = feasst.makePairLJCoulEwald(space, feasst.args(
            {"rCut" : str(rcut),
             "molType" : "none"}))
    elif input_parameters["pair_type"] == "LJ" or "WCA":
        pair = feasst.makePairLJ(space, feasst.args(
            {"rCut" : str(rcut),
             "molType" : "none"}))
    else:
        raise Exception('Unknown potential type: '+input_parameters["pair_type"])

    # open the json file created by process.py
    with open(in_file_name + "_out_types.json", 'r') as json_file:
        import json
        types = json.load(json_file)

    # initialize adsorbate, which must be first to use order parameter nMol0
    adsorbate_data_file = adsorbate_filename
    pair.initData(adsorbate_data_file)

    # initialize framework atom types
    # each atom type should have a corresponding data. LAMMPS-formatted file
    for atom in types:
        data_file_name = space.install_dir() + "/forcefield/data." + atom
        pair.initData(data_file_name)
        for i in range(types[atom]):
            pair.addMol(data_file_name)

    # read the framework coordinates sorted by process.py
    space.readXYZAlt(in_file_name + "_out.xyz")

    # initialize LJ cutoff
    pair.initAtomCut(1)
    pair.equateRcutForAllTypes()
    # initialize the tail correction
    if input_parameters["tail_type"] == "LFS":
        pair.linearShift(1) #Linear force shift tail
    elif input_parameters["pair_type"] == "WCA":
        pair.initWCA()
    # ewald Parameters
    if input_parameters["pair_type"] == "LJCoulEwald":
        alpha = input_parameters["Ewald"]["alpha"]
        k2max = input_parameters["Ewald"]["k2max"]
        pair.initKSpace(alpha, k2max)
        #pair.kmaxset(5,6,8)
        #pair.noErfTable() # turn off the error function table

    # initialize acceptance criteria
    criteria = feasst.makeCriteriaWLTMMC(feasst.args(
        {"beta" : str(1./(temp*feasst.idealGasConstant/1e3)),
         "activ" : str(activ),
         "mType" : "nmol0", # order parameter is the number of molecules of adsorbate
         "nMin" : str(0),
         "nMax" : str(1)}))
    for atom in types:
        criteria.addActivity(math.exp(-1))   # activity of framework doesn't matter, but error checked
    criteria.collectInit()
    criteria.tmmcInit()
    mc = feasst.WLTMMC(pair, criteria)

    # initialize total energy
    pair.initEnergy()
    peTot_bare = pair.peTot() #store the energy of the bare sorbent
    peLJ_bare = pair.peLJ()
    if input_parameters["pair_type"] == "LJCoulEwald":
        peQReal_bare = pair.peQReal()
        peQFrr_bare = pair.peQFrr()
        peQFrrSelf_bare = pair.peQFrrSelf()

    # add molecules up to nMolMin
    n_framework = space.nMol()
    mc.nMolSeek(n_framework + 1)
    n_atoms_tot = int(len(space.x())/space.dimen()) #is there a better way to do this?
    pair.initEnergy()
    # identify the "insertion" atoms
    mpart = [x for x in range(n_framework,n_atoms_tot)]

    # initialize the accumulators
    ncoeffs = input_parameters["ncoeffs"]
    Kcoeff = [ 0. ] * ncoeffs
    Kcoeff_var = [ 0. ] * ncoeffs
    factorial = [ np.math.factorial(j) for j in range(ncoeffs) ]
    count = 0
    rejected = 0

    #-----------------------------
    #  COEFFICIENT SCALING
    #   Note: Scale factors are simply to help avoid numerical overflow
    try:
        scale_coefficients = input_parameters["scale"]["active"]
    except:
        scale_coefficients = False
    if scale_coefficients:
        scale_factor = input_parameters["scale"]["factor"]**(-1.)
    else:
        scale_factor = 1.
    #-----------------------------

    # Adjustment term (bare sorbent Fourier-space terms)
    if input_parameters["pair_type"] == "LJCoulEwald":
        adjust = peQFrr_bare + (pair.peQFrrSelf()-peQFrrSelf_bare)
    else:
        adjust = 0.

    # Enable / Disable TQDM Progress Bar (Default: Enable Progress Bar)
    try:
        disable_flag = not input_parameters["progress_bar"]
    except:
        disable_flag = False

    # Monte Carlo Integration
    output_freq = int(input_parameters["trials"])/10
    for i in range(int(input_parameters["trials"])):
        space.randDisp(mpart,-1)
        space.randRotate(mpart,-1)

        deltaU = pair.multiPartEner(pyfeasst.list_to_int_vector(mpart), 1)-adjust

        count += 1
        for j in range(ncoeffs):
            if math.isnan(deltaU) or deltaU > 1.e19:
                value = 0.  #infinite energy
            else:
                value = np.exp( -criteria.beta() * deltaU ) * (deltaU*scale_factor)**j * ((-1.)**j)/factorial[j]
            # Accumulate scaled coefficients as a running average to avoid double precision overflow
            Kcoeff[j] += (value - Kcoeff[j])/float(count)
            Kcoeff_var[j] += (value**2 - Kcoeff_var[j])/float(count)
        if debug and (i+1) % output_freq == 0: print(i+1, deltaU, Kcoeff[0], Kcoeff[0]/float(count), rejected)
    # print()
    # Rescale coefficients & variance
    Kcoeff = [ Kcoeff[j]/(scale_factor**j) for j in range(ncoeffs) ]
    Kcoeff_var = [ (Kcoeff_var[j]/(scale_factor**(2*j)) - Kcoeff[j]**2) * float(input_parameters["trials"])/float(input_parameters["trials"]-1)
                   for j in range(ncoeffs) ]

    return Kcoeff, Kcoeff_var, criteria.beta()

# Wrapper function [preserves payload for the serial function]
def HenryCoefficient(debug=False,seed=123456,**kwargs):
    # Process the input XYZ file
    #  NOTE: this is done outside the main worker function so that the worker can
    #        be used by the parallel function.
    process(kwargs["adsorbent"])
    return HenryCoefficient_worker(debug=debug,seed=seed,**kwargs)

# Wrapper function [preserves payload for the serial function]
def HenryCoefficient_wrapper(arg):
    seed, kwargs = arg
    return HenryCoefficient_worker(seed=seed,**kwargs)

# Parallel Worker
def Parallel_HenryCoefficient(nthreads=4,seed=123456,**kwargs):
    input_dict = dict(kwargs.items())
    # Build a tuple of the arguments to pass to HenryCoefficient
    arg = [(seed+i, input_dict) for i in range(nthreads)]

    # Process the input XYZ file
    #  NOTE: this is done prior to the parallel fork to prevent collisions
    process(input_dict["adsorbent"])

    # Assemble and Execute the Parallel Job
    with mp.Pool(processes = nthreads) as pool:
        results = pool.map(HenryCoefficient_wrapper, arg)

    # Reassemble the parallel results
    ncoeffs = input_dict["ncoeffs"]
    trials_per_thread = input_dict["trials"]
    Kh = [0.] * ncoeffs
    Kh_var = [0.] * ncoeffs
    for coeffs, var, beta in results:
        for i in range(ncoeffs):
            Kh[i] += coeffs[i]*trials_per_thread
            Kh_var[i] += (var[i]*(trials_per_thread-1.)/trials_per_thread + coeffs[i]**2)*trials_per_thread

    Kh = [ Kh[i]/float(trials_per_thread*nthreads) for i in range(ncoeffs)]
    Kh_var = [ (Kh_var[i]/float(trials_per_thread*nthreads) - Kh[i]**2)
               *float(trials_per_thread*nthreads)/(float(trials_per_thread*nthreads)-1.)
                for i in range(ncoeffs)]

    return Kh, Kh_var, beta


if __name__ == "__main__":
    # execute only if run as a script
    temp_input = dict
    Kcoeffs, Kcoeffs_var, beta = HenryCoefficient(debug=True)
    for i,Kh in enumerate(Kcoeffs):
        print("Coefficient",i, ": ",Kh)
    print()
    print('Qst', 1./beta + Kcoeffs[1]/Kcoeffs[0]  )
