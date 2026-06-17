"""
In this tutorial, we simulate multiple temperatures and interpolate between them.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from feasst import fstio
from feasst import macrostate_distribution
from feasst import multistate_accumulator
import launch_04_lj_tm_parallel

def parse():
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/lj_new.txt',
          min_sweeps=10,
          cubic_side_length=8,
          window_alpha=2.5,
          hours_checkpoint=1,
          hours_terminate=5*24,
          num_jobs=32*3,
          collect_flatness=18,
          min_flatness=22)
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    params['script'] = __file__
    params['prefix'] = 'lj_ntrp'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['temps'] = [1.2, 1.1, 1]
    params['beta_mus'] = [-2.902929, -3.245568, -3.82307]
    # Choose max_particles relevant at the next simulated temp (for interpolation)
    params['nmaxs'] = [400, 425, 425]
    assert len(params['beta_mus']) == len(params['temps'])
    assert len(params['nmaxs']) == len(params['temps'])
    params['windows'] = list()
    assert params['num_jobs'] % len(params['temps']) == 0
    params['procs_per_temp'] = int(params['num_jobs']/len(params['temps']))
    for index,nm in enumerate(params['nmaxs']):
        params['windows'] += macrostate_distribution.window_exponential(
            alpha=params['window_alpha'], minimums=[0], maximum=nm,
            number=32, overlap=1, min_size=params['min_window_size'])
    return params, args

def sim_job_dependent_params(params):
    """ Define parameters that are dependent on the sim or job. """
    tsim = int(params['sim']/params['procs_per_temp'])
    params['beta'] = 1./params['temps'][tsim]
    params['beta_mu'] = params['beta_mus'][tsim]
    params['sim_start'] = params['job']*params['procs_per_job']
    params['sim_end'] = params['sim_start'] + params['procs_per_job'] - 1
    params['min_particles'] = params['windows'][params['sim']][0]
    params['max_particles'] = params['windows'][params['sim']][1]
    if params['max_particles'] < 125:
        params['min_sweeps'] = 50
    elif params['max_particles'] < 350:
        params['min_sweeps'] = 20

def post_process(params):
    import copy
    lnpis = list()
    V = params['cubic_side_length']**3
    for index,_ in enumerate(params['temps']):
    #for job in range(params['num_jobs']):
        #print('job', job)
        jstart = index*params['procs_per_temp']
        jstop = jstart + params['procs_per_temp'] - 1
        #print('jstart', jstart, 'jstop', jstop)
        lnpi = macrostate_distribution.splice_num_files(prefix=params['prefix']+'j', start=jstart, stop=jstop, suffix="s*_crit.csv", shift=False)
        #print('lnpi', lnpi.dataframe())
        en = multistate_accumulator.splice_num(prefix=params['prefix']+'j', start=jstart, stop=jstop, suffix='s*_en.csv')
        #print('en', en)
        lnpi.compute_single_component_derivatives(en)
        #print(lnpi.dataframe())
        lnpis.append(lnpi)
    temps = list(); rhov = list(); rhol = list()
    for tsim in range(len(params['temps']) - 1):
        #print('tsim', tsim)
        beta1 = 1./params['temps'][tsim]
        beta2 = 1./params['temps'][tsim+1]
        #print('temp1', 1/beta1, 'temp2', 1/beta2, 'beta1', beta1, 'beta2', beta2)
        #lnpis[tsim].plot(show=True)
        #lnpis[tsim+1].plot(show=True)
        lnpis[tsim].init_interpolate(lnpis[tsim+1], beta1=beta1, beta2=beta2)
        #print('A', lnpis[tsim]._ntrp_A)
        #print('xs', lnpis[tsim]._ntrp_xs)
        #betas = np.arange(beta1, beta2+1e-15, 0.01)
        num_ntrp = 100
        endpoint = False; num=num_ntrp
        if tsim == len(params['temps']) - 2:
            endpoint = True; num=num_ntrp + 1
        betas = 1/np.linspace(1./beta1, 1./beta2, num=num, endpoint=endpoint)
        #equil = pd.DataFrame(data={'temp': temps, 'rho_vap':None, 'rho_liq':None})
        for ibeta,beta in enumerate(betas):
            lnpin = copy.deepcopy(lnpis[tsim])
            #print('lnpin', lnpin.dataframe())
            lnpin.interpolate(beta)
            #print('lnpin', lnpin.dataframe()['ln_prob'])
            #lnpin.plot(show=True)
            lnpin.set_minimum_smoothing(30)
            lnpin.equilibrium()
            vap, liq = lnpin.split()
            temps.append(1/beta)
            rhov.append(vap.average_macrostate()/V)
            rhol.append(liq.average_macrostate()/V)
            print('T', temps[-1], 'rho_vap', rhov[-1], 'rho_liq', rhol[-1])
        #print(lnpis[tsim])
    equil = pd.DataFrame(data={"T":temps, "rho_vap":rhov, "rho_liq":rhol})
    equil = equil.iloc[::-1] # reverse
    srsw=pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/lj_srsw_equil.csv', comment='#')
    srsw=srsw[300:] # cut data down to the same temps
    #print(srsw)
    equil['rho_vap_srsw'] = srsw['rho_vap'].to_numpy()
    equil['rho_liq_srsw'] = srsw['rho_liq'].to_numpy()
    equil.to_csv(params['prefix']+'_equil.csv')
    #if True: # show plot
    if False: # do not show plot
        plt.scatter(equil['rho_vap'], equil['T'], color='red')
        plt.scatter(equil['rho_liq'], equil['T'], color='red', label='interpolated')
        plt.plot(equil['rho_vap_srsw'], equil['T'], color='black')
        plt.plot(equil['rho_liq_srsw'], equil['T'], color='black', label='SRSW')
        plt.xlabel(r'$\rho$', fontsize=16)
        plt.ylabel(r'$T$', fontsize=16)
        plt.legend(fontsize=16)
        plt.show()
    for phase in ['vap', 'liq']:
        diverged = equil[np.abs(equil['rho_'+phase] - equil['rho_'+phase+'_srsw']) > 0.02*equil['rho_'+phase+'_srsw']]
        if len(diverged) != 0:
            print('diverged', diverged)
            assert False

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_job_dependent_params=sim_job_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_job,
                          args=arguments)
