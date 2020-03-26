# Import libraries
from lib.time.time import dt
from numpy import sqrt, sin, cos, arctan2

# Satellite ECEF position
def satpos(ttr,toe,ROOTa,DELTAn,M0,e,omega,Cus,Cuc,Crs,Crc,Cis,Cic,i0,iDOT,OMEGA0,OMEGADOT):

    # IS-GPS-200K constants(WGS84)
    pi = 3.1415926535898            # WGS84 version of pi
    c = 2.99792458e8                # Speed of Light [m/s]
    OMEGADOTe = 7.2921151467e-5     # Earth's rotation rate [rad/s]
    GM = 3.986005e14                # Earth's gravitational constant [m^3/s^2]

    # Anomalies of the Keplerian orbit
    a = ROOTa**2                    # Semi-major axis [m]
    n0 = sqrt(GM/a**3)              # Mean angular velocity [rad/sec]
    t = dt(ttr, toe)                # Time from reference epoch [s]
    n = n0 + DELTAn                 # Corrected mean motion [rad/s]
    M = M0 + n*t                    # Mean anomaly [rad]

    # Eccentric anomaly
    epsilon = 1e-10
    E = M, E0 = 0

    while abs(E - E0) > epsilon:
        E0 = E
        E = M + e*sin(E0)

    # True anomaly
    v = arctan2(sqrt(1 - e**2)*sin(E), cos(E) - e)
