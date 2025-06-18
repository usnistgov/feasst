"""
In this tutorial, we simulate multiple temperatures and interpolate between them.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution
from pyfeasst import multistate_accumulator
import launch_04_lj_tm_parallel

def parse():
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/lj_new.txt',
          min_sweeps=5,
          cubic_side_length=8,
          window_alpha=2.5,
          hours_checkpoint=1,
          hours_terminate=5*24,
          num_nodes=3,
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
    assert len(params['nmaxs']) == params['num_nodes']
    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    for node in range(params['num_nodes']):
        if node == params['node']:
            params['beta'] = 1./params['temps'][node]
            params['beta_mu'] = params['beta_mus'][node]
            params['max_particles'] = params['nmaxs'][node]
    params['sim_start'] = params['node']*params['procs_per_node']
    params['sim_end'] = params['sim_start'] + params['procs_per_node'] - 1
    params['windows'] = macrostate_distribution.window_exponential(
        alpha=params['window_alpha'], minimums=[params['min_particles']], maximum=params['max_particles'],
        number=params['procs_per_node'], overlap=1, min_size=params['min_window_size'])
    params['min_particles'] = params['windows'][params['sim']-params['sim_start']][0]
    params['max_particles'] = params['windows'][params['sim']-params['sim_start']][1]

def post_process(params):
    import copy
    lnpis = list()
    V = params['cubic_side_length']**3
    for node in range(params['num_nodes']):
        #print('node', node)
        lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n'+str(node)+'s', suffix='_crit.csv', shift=False)
        #print('lnpi', lnpi.dataframe())
        en = multistate_accumulator.splice(prefix=params['prefix']+'n'+str(node)+'s', suffix='_en.csv')
        #print('en', en)
        lnpi.compute_single_component_derivatives(en)
        #print(lnpi.dataframe())
        lnpis.append(lnpi)
    temps = list(); rhov = list(); rhol = list()
    for node in range(params['num_nodes'] - 1):
        #print('node', node)
        beta1 = 1./params['temps'][node]
        beta2 = 1./params['temps'][node+1]
        #print('beta1', beta1, 'beta2', beta2)
        #lnpis[node].plot(show=True)
        lnpis[node].init_interpolate(lnpis[node+1], beta1=beta1, beta2=beta2)
        #betas = np.arange(beta1, beta2+1e-15, 0.01)
        num_ntrp = 100
        endpoint = False; num=num_ntrp
        if node == params['num_nodes'] - 2: endpoint = True; num=num_ntrp + 1
        betas = 1/np.linspace(1./beta1, 1./beta2, num=num, endpoint=endpoint)
        #equil = pd.DataFrame(data={'temp': temps, 'rho_vap':None, 'rho_liq':None})
        for ibeta,beta in enumerate(betas):
            lnpin = copy.deepcopy(lnpis[node])
            lnpin.interpolate(beta)
            #lnpin.plot(show=True)
            lnpin.set_minimum_smoothing(30)
            lnpin.equilibrium()
            vap, liq = lnpin.split()
            temps.append(1/beta)
            rhov.append(vap.average_macrostate()/V)
            rhol.append(liq.average_macrostate()/V)
            print('T', temps[-1], 'rho_vap', rhov[-1], 'rho_liq', rhol[-1])
        #print(lnpis[node])
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
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
