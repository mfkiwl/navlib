# Import libraries
from numpy import array, eye, vstack, sin, cos, sqrt, arctan2
from numpy.linalg import norm


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


# Rotation matrix to Euler angles
# R = Rz(yaw)@Ry(pitch)@R(roll)
# Note: singular for pitch = 0 (Gimbal lock)
# DCM2Euler
def DCM2Euler(C):
    roll = arctan2(C[2, 1], C[2, 2])
    pitch = arctan2(-C[2, 0], sqrt(C[2, 1]**2 + C[2, 2]**2))
    yaw = arctan2(C[1, 0], C[0, 0])
    return roll, pitch, yaw


# Euler2DCM
def Euler2DCM(roll, pitch, yaw):
    return Rz(yaw)@Ry(pitch)@Rx(roll)


# Skew matrix
def skew(x):
    return array([[0, -x[2, 0], x[1, 0]],
                  [x[2, 0], 0, -x[0, 0]],
                  [-x[1, 0], x[0, 0], 0]])


# Axis-angle to DCM (Rodriguez formula)
def axis_ang2DCM(theta, r):
    Sr = skew(r)
    return eye(3) + sin(theta)*Sr + (1 - cos(theta))*Sr@Sr


# Quaternion multiplication
def qmult(p, q):
    return array([[p[0, 0]*q[0, 0] - p[1, 0]*q[1, 0] - p[2, 0]*q[2, 0] - p[3, 0]*q[3, 0]],
                  [p[1, 0]*q[0, 0] + p[0, 0]*q[1, 0] - p[3, 0]*q[2, 0] + p[2, 0]*q[3, 0]],
                  [p[2, 0]*q[0, 0] + p[3, 0]*q[1, 0] + p[0, 0]*q[2, 0] - p[1, 0]*q[3, 0]],
                  [p[3, 0]*q[0, 0] - p[2, 0]*q[1, 0] + p[1, 0]*q[2, 0] + p[0, 0]*q[3, 0]]])


# Quaternion conjugate
def qconj(q):
    return vstack([q[0], -q[1:]])


# Quaternion inverse
def qinv(q):
    return qconj(q)/norm(q)**2


# DCM to quaternion
def DCM2quat(C):
    q = 1/2*sqrt(1 + C[0, 0] + C[1, 1] + C[2, 2])
    return array([[q],
                  [-1/(4*q)*(C[1, 2] - C[2, 1])],
                  [-1/(4*q)*(C[2, 0] - C[0, 2])],
                  [ 1/(4*q)*(C[0, 1] - C[1, 0])]])


# quaternion to DCM
def quat2DCM(q):
    return array([[ q[0, 0]**2 + q[1, 0]**2 - q[2, 0]**2 - q[3, 0]**2, 2*(q[1, 0]*q[2, 0] + q[3, 0]*q[0, 0]), -2*(q[1, 0]*q[3, 0] - q[2, 0]*q[0, 0])],
                  [ 2*(q[1, 0]*q[2, 0] - q[3, 0]*q[0, 0]), q[0, 0]**2 + q[2, 0]**2 - q[1, 0]**2 - q[3, 0]**2, -2*(q[2, 0]*q[3, 0] + q[0, 0]*q[1, 0])],
                  [-2*(q[1, 0]*q[3, 0] + q[2, 0]*q[0, 0]), -2*(q[2, 0]*q[3, 0] - q[0, 0]*q[1, 0]), q[0, 0]**2 + q[3, 0]**2 - q[1, 0]**2 - q[2, 0]**2]])


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
    roll = arctan2(ay/az)
    pitch = arctan2(ax/az)
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
