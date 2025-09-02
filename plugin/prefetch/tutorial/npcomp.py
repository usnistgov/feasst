import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear_fit(x, b):
    return -0.5*x + b

def post_process(params):
    z_factor = 3
    sims={'dr':['np1/', 'np2/'], 'np':[1, 2]}
    for index,dr in enumerate(sims['dr']):
        time = pd.read_csv(dr+'ljgibbs000_cpu.csv', delim_whitespace=True, header=None)
        print(time[1])
        en = pd.read_csv(dr+'ljgibbs000_liquid_en.csv', header=None, comment="a")
        print(en[0])
        equil=50
        logt = np.log(time[1][equil:])
        print('logt', logt)
        logs = np.log(en[2][equil:])
        print('logs', logs)
        popt, pcov = curve_fit(linear_fit, logt, logs)
        plt.scatter(np.log(time[1]), np.log(en[2]))
        plt.plot(logt, linear_fit(logt, popt[0]), label=dr+" "+str(popt[0]))
        print(popt[0])
        #plt.xscale('log')
        #plt.yscale('log')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    post_process(None)
