import numpy as np
import math

__author__ = 'Pierre Winter'

"""
This module calculates sums and densities of vibrational states for a molecule at a given energy E.
Counting algorithms with various levels of approximations are implemented as functions:

cl_count - Classical approximation, only useful at high energies
MR_count - Semi-classical Marcus-Rice approximation, considers zero-point energy
WR_count - Semi-classical Whitten-Rabinovitch approximation, considers a fraction of the zero-point energy
BS_count - Beyer-Swinehart exact counting algorithm, gold standard in chemistry
"""


class state_count(object):
    def __init__(self, E, freqs, TS=False, linear=False, ref='zpe'):
        """
        This module calculates sums and densities of states for a molecule at a given energy E.
        Currently only sums and densities of vibrational states are implemented.

        :param E: Total energy in the system
        :param freqs: Array of frequencies in wavenumbers (cm-1)
        :param TS: Flag which determines if molecule is a minimum or a TS geometry
        :param linear: Flag which indicates if a molecule has a linear geometry
        :param ref: Defines the reference of energy, can be from bottom of potential well ('bot')
        or can be from the zero-point energy of the molecule ('zpe')
        """

        self.E = E  # units of cm-1
        self.freqs = np.sort(freqs)  # units of cm-1
        self.hbar = 1.  # set to 1 since Planck's constant is already taken into account inside frequencies
        self.ref = ref
        if TS:
            assert self.freqs[0] < 0
            self.freqs = self.freqs[1:]
        if linear:
            self.vibs = np.array(self.freqs[5:])
            self.sum = np.sum(self.vibs)
            self.n_osc = len(self.vibs)
        else:
            self.vibs = np.array(self.freqs[6:])
            self.sum = np.sum(self.vibs)
            self.n_osc = len(self.vibs)

    def zpe(self):
        return self.hbar * self.sum / 2

    def cl_count(self):
        """ Calculates a classical sum and density of states using statistical mechanics. """

        if self.ref == 'zpe':
            E = self.E
        else:
            E = self.E - self.zpe()

        s = self.n_osc
        prod = np.prod(self.vibs)  # No need to multiply vibs by Planck's constant h because vibs in cm-1 already
        summ = E**s / (prod * math.factorial(s))
        dens = summ * s / E

        return summ, dens

    def MR_count(self):
        """ Calculates a semi-classical sum of states using Marcus-Rice approach. """

        if self.ref == 'zpe':
            E = self.E
        else:
            E = self.E - self.zpe()

        s = self.n_osc
        prod = np.prod(self.vibs)  # No need to multiply vibs by Planck's constant h because vibs in cm-1 already
        summ = (E + self.zpe())**s / (prod * math.factorial(s))
        dens = summ * (s / (E + self.zpe()))

        return summ, dens

    def WR_count(self):
        """ Calculates a semi-classical sum of states using Whitten-Rabinovitch approach. """

        if self.ref == 'zpe':
            E = self.E
        else:
            E = self.E - self.zpe()

        s = self.n_osc
        denom = np.prod(self.vibs) * math.factorial(s)
        vratio = (np.sum(self.vibs**2) / s) / (self.sum / s)**2
        b = vratio * (s - 1) / s
        eprime = E / self.zpe()
        if 0.01 < eprime < 1.0:
            w = 1. / (5.00 * eprime + 2.73 * eprime**0.50 + 3.51)
            dwde = -1. / (5.00 * eprime + 2.73 * eprime**0.50 + 3.51)**2 * (5 + 1.365*eprime**(-0.5))
        elif 1.0 <= eprime:
            w = np.exp(-2.4191 * eprime**0.25)
            dwde = np.exp(-2.4191 * eprime**0.25) * (-0.604775 * eprime**(-0.75))
        else:  # eprime < 0.1
            print('\nWARNING: value of eprime is too small at %g\n' % eprime)
            w = 0
            dwde = 0

        a = 1 - b * w
        num = (E + a * self.zpe())**s
        summ = num / denom

        deriv = 1 - b * dwde
        dens = summ * (s / (E + a * self.zpe()) * deriv)

        return summ, dens

    def BS_count(self, Egrain):
        """ Calculates an exact sum of states using the general Bayer-Swinehart algorithm.
        Egrain defines the spacing within which the counts are taken in units of wavenumbers (cm-1).
        Egrain=50 is usually a safe spacing for molecular vibrations."""

        if self.ref == 'zpe':
            E = self.E
        else:
            E = self.E - self.zpe()

        R = [int(round(v/Egrain)) for v in self.vibs]
        T = np.zeros(100000)
        M = len(T)-1
        T[0] = 1
        for r in R:
            T[r] += T[0]
            T[r+1] += T[1]
            for i in range(r+1, M):
                T[i] += T[i-r]
            T[M] += T[M-r]
        RM = int(round(E / Egrain))
        summ = int(np.sum(T[:RM]))
        dens = T[RM-1] / Egrain

        return summ, dens


if __name__ == "__main__":
    """
    This is an example of how to use the module.
    We calculate the sums and densities of states for the cylcopropane molecule.
    """

    energy = 52410  # in units of wavenumbers (cm-1)
    freqs_cyclopropane = np.array([3221, 3221, 3221, 3221, 3221, 3221,
                                 1478, 1478, 1478,
                                 1118, 1118, 1118, 1118, 1118, 1118, 1118,
                                 878, 878, 878,
                                 749, 749,
                                 0, 0, 0, 0, 0, 0], dtype=float)  # in units of wavenumbers (cm-1)
    count_vibs = state_count(E=energy, freqs=freqs_cyclopropane, TS=False, linear=False, ref='zpe')

    cl_sum, cl_dens = count_vibs.cl_count()
    mr_sum, mr_dens = count_vibs.MR_count()
    wr_sum, wr_dens = count_vibs.WR_count()
    bs_sum, bs_dens = count_vibs.BS_count(Egrain=50)  # Egrain in units of wavenumbers (cm-1)

    print('classical\n  sum: {:.3e},  dens: {:.3e}'.format(cl_sum, cl_dens))
    print('marcus-rice\n  sum: {:.3e},  dens: {:.3e}'.format(mr_sum, mr_dens))
    print('whitten-rabinovitch\n  sum: {:.3e},  dens: {:.3e}'.format(wr_sum, wr_dens))
    print('beyer-swinehart\n  sum: {:.3e},  dens: {:.3e}'.format(bs_sum, bs_dens))

