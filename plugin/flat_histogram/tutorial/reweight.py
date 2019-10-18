import math
import copy
import numpy as np

class reweight():
    """Basic histogram reweighting functionality can be used to obtain information from the
    macrostate distribution
    """
    def __init__(self):
        self.set_num_smooth()
        self.set_phase_boundary()
        self.set_lnzs_guess(-2)
        self._saturation_found = False
        self._lnpi_vap_sat = None
        self._lnpi_liq_sat = None
        self.set_lnpi(None)
        self.set_beta()

    def set_num_smooth(self, num_smooth=20):
        """Find global minimum from this many macrostates"""
        self._num_smooth = num_smooth

    def set_phase_boundary(self, phase_boundary=-1):
        """Set the phase boundary"""
        if phase_boundary == -1:
            self._automatic_phase_boundary = True
            self._phase_boundary = -1
        else:
            self._automatic_phase_boundary = False
            self._phase_boundary = phase_boundary

    def get_phase_boundary(self):
        """Get the last phase boundary"""
        return self._phase_boundary

    def set_lnzs_guess(self, lnzsguess):
        """Set the guess for the activity(lnz) of saturation."""
        self._lnzsguess = lnzsguess

    def set_lnpi(self, lnpi):
        """Set the natural logarithm of the macrostate distribution"""
        self._lnpi = lnpi

    def set_beta(self, beta=1./1.2):
        """Set the inverse temperature"""
        self._beta = beta

    def _local_min_smooth(self, vec):
        """Obtain local minima with smoothing"""
        mins = []
        for index in range(1, len(vec)-1):
            lower = index-self._num_smooth
            if lower < 0:
                lower = 0
            upper = index+self._num_smooth
            if upper >= len(vec):
                upper = len(vec) - 1
            minimum = np.amin(vec[lower:upper+1])
            if math.fabs(minimum - vec[index]) < 1e-7:
                mins.append(index)
        return mins

    def _nboundary(self, lnpi):
        """Return the number of particles which represents the phase boundary"""
        res = self._local_min_smooth(lnpi)
        if res:
            return res[0]
        return None

    def area(self, lnpi):
        """Return the area"""
        sum_area = 0
        for _, lnpii in enumerate(lnpi):
            sum_area += math.exp(lnpii)
        return sum_area

    def norm(self, lnpi):
        """Return the normalized log Pi"""
        lnpi = np.array(lnpi)
        lnpi_working = lnpi
        shift = math.log(self.area(lnpi_working))
        lnpi -= shift
        return lnpi

    def reweight(self, lnpi, delta_mu):
        for i, _ in enumerate(lnpi):
            lnpi[i] += i*self._beta*delta_mu

    def _sqprobdiff(self, mu, murw):
        """Separate the phase boundary and return the squared difference in the areas of these
        macrostate distributions
        """
        lnpi = copy.deepcopy(self._lnpi)
        self.reweight(lnpi, murw - mu)
        if self._automatic_phase_boundary:
            boundary = self._nboundary(lnpi)
        else:
            boundary = self._phase_boundary
        # if no phase boundary, encourge multiple peaks by biasing toward similar
        # orders of magnitude of the first and last value of the lnpi
        if boundary is None:
            return (lnpi[0] - lnpi[-1])**2
        vap = lnpi[:boundary]
        liq = lnpi[boundary:]
        self._phase_boundary = boundary
        self._lnpi_vap_sat = vap
        self._lnpi_liq_sat = liq
        return (math.log(self.area(vap)) - math.log(self.area(liq)))**2

    def lnzsat(self, lnpi, mu):
        """Find saturation by equating the probabilities of the two phases"""
        from scipy.optimize import minimize
        lnzsguess = self._lnzsguess
        self._lnpi = copy.deepcopy(lnpi)
        res = minimize(lambda lnzrwi: self._sqprobdiff(mu, murw=lnzrwi[0]), lnzsguess, tol=1e-8)
        self._saturation_found = True
        return res["x"][-1]

    def average_macrostate(self, lnpi, first_macro=0):
        """Return the average macrostate"""
        average = 0.
        for index, value in enumerate(lnpi):
            average += (index + first_macro)*math.exp(value)
        return average/self.area(lnpi)

    def saturation_properties(self, volume, lnpi, mus):
        """Return saturation properties"""
        assert self._saturation_found
        return_dict = {"mu_sat": mus}
        return_dict["phase_boundary"] = self.get_phase_boundary()
        return_dict["vapor"] = {"pressure":
            (-lnpi[0] + math.log(self.area(self._lnpi_vap_sat)))/volume/self._beta,
            "density": self.average_macrostate(self._lnpi_vap_sat)/volume}
        return_dict["liquid"] = {"pressure":
            (-lnpi[0] + math.log(self.area(self._lnpi_liq_sat)))/volume/self._beta,
            "density": self.average_macrostate(self._lnpi_liq_sat, len(self._lnpi_vap_sat))/volume}
        return return_dict
