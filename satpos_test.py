# Import libraries
from numpy import array


# Receiver position [m]
Xr = array([[3172870.7170],
            [604208.2810],
            [5481574.2300]])

# Satellite G01 broadcast ephemerides (RINEX)
ttr = 8134                        # [s]
toe = 7200                        # [s]
ROOTa = 5.153634706497e+03        # [sqrt(m)]
DELTAn = 4.646979279625e-09       # [rad/s]
M0 = 9.760178388778e-01           # [rad]
e = 9.364774916321e-03            # [unitless]
omega = 7.546597134633e-01        # [rad]
Cus = 1.266598701477e-07          # [rad]
Cuc = -2.680346369743e-06         # [rad]
Crs = -5.456250000000e+01         # [m]
Crc = 3.865625000000e+02          # [m]
Cis = 1.285225152969e-07          # [rad]
Cic = -8.940696716309e-08         # [rad]
i0 = 9.785394956406e-01           # [rad]
iDOT = -3.500145795122e-10        # [rad/s]
OMEGA0 = -1.328259931335e+00      # [rad]
OMEGADOT = -8.668218208939e-09    # [rad/s]
