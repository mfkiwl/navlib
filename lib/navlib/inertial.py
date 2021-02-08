# Import libraries
from numpy import array
from scipy import sin, cos, sqrt, arctan


# Quaternion
def DCM2quat(C):
    q = 1/2*sqrt(1 + C[0, 0] + C[1, 1] + C[2, 2])

    return array([[q],
                  [1/(4*q)*(C[1, 2] - C[2, 1])],
                  [1/(4*q)*(C[2, 0] - C[0, 2])],
                  [1/(4*q)*(C[0, 1] - C[1, 0])]])


def quat2DCM(q):

    return array([[q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 2*(q[1]*q[2] + q[3]*q[0]), 2*(q[1]*q[3] - q[2]*q[0])],
                  [2*(q[1]*q[2] - q[3]*q[0]), q[0]**2 + q[2]**2 - q[1]**2 - q[3]**2, 2*(q[2]*q[3] + q[0]*q[1])],
                  [2*(q[1]*q[3] + q[2]*q[0]), 2*(q[2]*q[3] - q[0]*q[1]), q[0]**2 + q[3]**2 - q[1]**2 - q[2]**2]])


# Directional Cosine Matrix (DCM)
# x-axis rotation
def Rx(rx):

    return array([[1, 0, 0],
                  [0, cos(rx), -sin(rx)],
                  [0, sin(rx),  cos(rx)]])


# y-axis rotation
def Ry(ry):

    return array([[cos(ry), 0, sin(ry)],
                  [0, 1, 0],
                  [-sin(ry), 0, cos(ry)]])


# z-axis rotation
def Rz(rz):

    return array([[cos(rz), -sin(rz), 0],
                  [sin(rz),  cos(rz), 0],
                  [0,        0,       1]])


# Rotate from e-frame to g-frame (ned)
def Ce_g(lat, lon):

    return array([[-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
                  [-sin(lon), cos(lon), 0],
                  [-cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)]])


# Rotate from b-frame to g-frame
def Cb_g(roll, pitch, yaw):

    return Rz(yaw)@Ry(pitch)@Rx(roll)


# Estimate roll and pitch from acceleration (ned)
def align(ax, ay, az):
    roll = arctan(ay/az)
    pitch = arctan(ax/az)

    return roll, pitch


# Coordinate axis definitions
ned2enu = array([[0, 1, 0],
                 [1, 0, 0],
                 [0, 0, -1]])

nwu2enu = array([[0, -1, 0],
                 [1, 0, 0],
                 [0, 0, 1]])

nwu2ned = array([[1, 0, 0],
                 [0, -1, 0],
                 [0, 0, -1]])
