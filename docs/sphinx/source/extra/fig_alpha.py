import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np


def alpha(T):
    # Values for cpx in MinDB.csv
    alpha0 = 0.000053
    alpha1 = 5.92E-09
    alpha2 = -0.0122
    alpha3 = 0.672
    return alpha0 + alpha1*T + alpha2/T + alpha3/T/T


d = np.loadtxt(
    '../../../../VelocityConversion/AlphaDB.csv',
    delimiter=';',
    skiprows=3,
    usecols=[0, 1, 2]
)
tri = Triangulation(d[:, 1], d[:, 0])

T_min = np.unique(d[:, 1].min())
T_max = 2000
P_max = 8

alpha_fmt = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*1e6))

fig, ax = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'wspace': 0.5})

ax[0].plot(np.unique(d[:, 1]), alpha(np.unique(d[:, 1])), label=r'$\alpha(T)$')
ax[0].plot([T_min, T_max], [0.000053]*2, label=r'$\alpha =$ const.')
ax[0].set_xlim(T_min, T_max)
ax[0].yaxis.set_major_formatter(alpha_fmt)
ax[0].set_xlabel('Temperature / K')
ax[0].set_ylabel('Isothermal expansion coefficient / $10^{-6}$ K$^{-1}$')
ax[0].set_title(r'$\alpha = $const. and $\alpha(T)$')
ax[0].legend()

c = ax[1].tricontourf(tri, d[:, 2]*1e6)
ax[1].set_xlim(T_min, T_max)
ax[1].set_xlabel('Temperature / K')
ax[1].set_ylabel('Pressure / GPa')
ax[1].set_ylim(top=P_max)
ax[1].set_title(r'$\alpha(P,T)$')

fig.colorbar(c, ax=ax[:], orientation='vertical',
             label='Isothermal expansion coefficient / $10^{-6}$ K$^{-1}$')

fig.savefig('alpha.png', bbox_inches='tight')
