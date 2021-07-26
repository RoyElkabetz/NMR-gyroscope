import numpy as np

# PHYSICAL CONSTANTS
C = 2.99792458 * 1e8        # [m / s]

# SPINS
I85 = 5 / 2                 # Rubidium 85 Nuclear Spin (I)
I87 = 3 / 2                 # Rubidium 87 Nuclear Spin (I)

# GYROMAGNETIC RATIOS
G129 = -74.4069 * 1e6       # Xenon129 [rad / s / T]  http://www.acadiau.ca/~bellis/resources/nmr/isotopes.html
G131 = 22.0564 * 1e6        # Xenon131 [rad / s / T]  http://www.acadiau.ca/~bellis/resources/nmr/isotopes.html
Ge = 1.760859644 * 1e11     # isolated electron [rad / s / T]  https://en.wikipedia.org/wiki/Gyromagnetic_ratio
G85 = Ge / (2 * I85 + 1)    # Rubidium 85 ESR gyromagnetic ratio [rad / s / T]
G87 = Ge / (2 * I87 + 1)    # Rubidium 85 ESR gyromagnetic ratio [rad / s / T]

# UNITS CONVERSIONS
G2T = 1e-4                  # Gauss to Tesla conversion. i.e. 5 Gauss * G2T = 5e-4 Tesla
T2G = 1e4                   # Tesla to Gauss conversion. i.e. 5 Tesla * T2G = 5e4 Gauss
D2R = np.pi / 180           # degrees to radians conversion. i.e. 90 deg * D2R = 1.57 [rad]
R2D = 180 / np.pi           # radians to degrees conversion. i.e. 2pi * R2D = 360 [deg]

# SYMBOLS
PI = chr(960)
PHI = chr(966)
OMEGA = chr(969)
BIG_OMEGA = chr(937)
