"""
This module computes quantities to compare with scattering experiments.
"""

import math
import numpy as np
import pandas as pd

def structure_factor(distances, radial_distribution, frequency, number_density):
    """
    Return the structure factor at a given frequency by Fourier transform.
    Assume that the distances correspond with the radial distribution and are evenly spaced.
    Each distance is the center of the bin, and the width is twice the first distance.

    :param float number_density: the number of particles divided by the volume.

    >>> from pyfeasst import scattering
    >>> gr = pd.read_csv('../../tests/gr.csv')
    >>> sq = scattering.structure_factor(
    ...     distances=gr['r'],
    ...     radial_distribution=gr['g_HH'],
    ...     frequency=0.4,
    ...     number_density=530/90**3)
    >>> round(sq, 8)
    0.59292662
    """
    struct_fac = 0.
    assert len(distances) == len(radial_distribution)
    bin_width = 2*distances[0]
    for index, distance in enumerate(distances):
        qr = distance*frequency
        struct_fac += bin_width*distance*distance*math.sin(qr)*(radial_distribution[index] - 1)/qr
    return 1 + 4*math.pi*number_density*struct_fac

def gen_ff_file(qs, sigmas, output_file):
    """
    Generate a file with the form factors of spheres with the given sigma(diameter) values.

    >>> from pyfeasst import scattering
    >>> qs = [0.447022, 0.48368, 0.0698132]
    >>> scattering.gen_ff_file(qs=qs,
    ...                        sigmas=[1.5200, 5.3283, 4.5405, 4.4130, 4.6319],
    ...                        output_file='ff.csv')
    >>> import pandas as pd
    >>> df = pd.read_csv('ff.csv')
    >>> round(df['q'][0], 6) == qs[0]
    True
    >>> round(df['p0'][0], 4) == 14.0422
    True
    >>> round(df['p0'][1], 4) == 13.9303
    True
    >>> round(df['p4'][0], 3) == 263.059
    True
    >>> round(df['p4'][1], 3) == 241.494
    True
    >>> round(df['p0'][2], 4) == 14.6937
    True
    >>> round(df['p4'][2], 3) == 411.925
    True
    """

    df = pd.DataFrame(qs, columns=['q'])
    for isigma, sigma in enumerate(sigmas):
        ff = np.zeros(len(qs))
        for iq, q in enumerate(qs):
            volume = 4*np.pi*np.power(sigma, 3)/3.
            ff[iq] = volume*3*(np.sin(q*sigma)-q*sigma*np.cos(q*sigma))/np.power(q*sigma, 3)
        df['p'+str(isigma)] = ff
    df.to_csv(output_file)

def intensity(gr_file, num_density, iq_file=None, pqs=None, skip=1, rdf_shift=0.05):
    """
    Return the scattering intensity given a radial distribution function file with all rdfs
    between all site types, as well as intramolecular histograms.

    Use an iq_file output by feasst::Scattering to obtain the frequencies and scattering lengths.
    Fourier transform the grs to obtain the intensity.

    :param Pandas DataFrame pqs: instead of iq_file, provide pq directly as dataframe.
    :param int skip: compute this many less frequencies than are present in iq_file.
    :param float rdf_shift: this percent of the highest distance are shifted to one.

    >>> from pyfeasst import scattering
    >>> df = intensity('../../tests/grn30.csv', iq_file='../../tests/iq.csv', num_density=3/90**3, skip=100)
    >>> round(df['q'][0], 5)
    0.06981
    >>> round(df['iq'][0], 8)
    73049.60991857
    """
    rdf = pd.read_csv(gr_file, comment="#")
    dr = 2.*rdf['r'][0]

    # obtain the number of site types
    file1 = open(gr_file, 'r')
    lines = file1.readlines()
    file1.close()
    exec('iprm={' + lines[0][1:] + '}', globals())
    num_site_types = iprm['num_site_types']

    # shift
    if rdf_shift > 0 and rdf_shift < 1:
      for itype in range(num_site_types):
        for jtype in range(itype, num_site_types):
          gij = 'g' + str(itype) + '-' + str(jtype)
          cut = int((1-rdf_shift)*len(rdf[gij]))
          rdf[gij] /= np.average(rdf[gij][cut:])

    if pqs is None:
        assert iq_file is not None
        iq = pd.read_csv(iq_file, comment="#")
        iqgrp = iq.groupby('q', as_index=False)
        iq = iqgrp.mean()
    else:
        assert iq_file is None
        iq = pqs

    iq_intra = np.zeros(1+int((len(iq['q'])-1)/skip))
    iq_inter = np.zeros(len(iq_intra))
    iq_ff = np.zeros(len(iq_intra))
    iq_q = np.zeros(len(iq_intra))
    for iq_index, q in enumerate(iq['q']):
      if iq_index % skip == 0:
        #print('progress %:', 100*iq_index/len(iq['q']))
        skip_index = int(iq_index/skip)

        # form factor
        iq_ff_acc = 0.
        for itype in range(num_site_types):
          iq_ff_acc += iq['p'+str(itype)].values[iq_index]**2
        iq_ff[skip_index] = iq_ff_acc

        # intra and inter
        iq_intra_acc = 0.
        iq_inter_acc = 0.
        for gr_index, r in enumerate(rdf['r']):
          sqrinvqr = math.sin(q*r)/q/r
          rmin = dr*gr_index
          rmax = dr*(gr_index + 1)
          dv = 4./3.*math.pi*(rmax**3 - rmin**3)
          ideal_num = dv*num_density
          for itype in range(num_site_types):
            for jtype in range(itype, num_site_types):
              pqi = iq['p'+str(itype)].values[iq_index]
              pqj = iq['p'+str(jtype)].values[iq_index]
              factor = pqi*pqj*sqrinvqr
              if itype == jtype:
                factor *= 2
              iq_intra_acc += factor*rdf['h'+str(itype)+'-'+str(jtype)].values[gr_index]
              iq_inter_acc += factor*(rdf['g'+str(itype)+'-'+str(jtype)].values[gr_index] - 1)*ideal_num
        iq_q[skip_index] = q
        iq_intra[skip_index] = iq_intra_acc
        iq_inter[skip_index] = iq_inter_acc

    for iq in [iq_ff, iq_intra, iq_inter]:
        iq /= num_site_types**2
    return pd.DataFrame({'q': iq_q, 'iq': iq_ff + iq_intra + iq_inter, 'ff': iq_ff, 'intra': iq_intra, 'inter': iq_inter})

if __name__ == "__main__":
#    for num in [30]:
#    #for num in [3, 30]:
#      df = intensity('../../../plugin/chain/tutorial/test/cg7mab_gr_n'+str(num)+'.csv', '../../../plugin/chain/tutorial/test/cg7mab_iq_n'+str(num)+'.csv', num_density=num/90**3, skip=1)
#      print(df)
#      df.to_csv('tt'+str(num)+'.csv')
#    df = intensity('../../../plugin/chain/tutorial/cg7mab_gr_n30.csv', '../../../plugin/chain/tutorial/cg7mab_iq_n30.csv', num_density=30/90**3, skip=10)
    #df = intensity('../../tests/grn30.csv', '../../tests/iqn30.csv', num_density=30/90**3)
    import doctest
    doctest.testmod()
