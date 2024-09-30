# script_name.py
#
# Usage: python3 script_name.py arg1 arg2 ...
# Text explaining script usage
# Parameters:
# arg1: description of argument 1
# arg2: description of argument 2
# ...
# Output:
# A description of the script output
#
# Written by Erik Judy
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.
# import Python modules
# e.g., import math # math module
import sys
import math
import numpy as np
# "constants"
# e.g., R_E_KM = 6378.137
R_E_KM = 6378.137
E_E = 0.081819221456
# helper functions
## function description
o_x_km=float('nan')
o_y_km=float('nan')
o_z_km=float('nan')
x_km=float('nan')
y_km=float('nan')
z_km=float('nan')
def calc_denom(ecc,lat_rad):
    return math.sqrt(1.0-ecc**2.0*math.sin(lat_rad)**2.0)

# pass
# initialize script arguments
# arg1 = '' # description of argument 1
# arg2 = '' # description of argument 2
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km  = float(sys.argv[3])
    x_km = float(sys.argv[4])
    y_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print(\
     'Usage: '\
     'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
    )
    exit()
# exit()
# write script below this line
ecef_vector=np.array([o_x_km-x_km, o_y_km-y_km, o_z_km-z_km])
lon_rad = math.atan2(y_km,x_km)
lon_deg = lon_rad*180.0/math.pi
# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(z_km/math.sqrt(x_km**2+y_km**2+z_km**2))
r_lon_km = math.sqrt(x_km**2+y_km**2)
prev_lat_rad = float('nan')
# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
    denom = calc_denom(E_E,lat_rad)
    c_E = R_E_KM/denom
    prev_lat_rad = lat_rad
    lat_rad = math.atan((z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
    count = count+1
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E
lat_deg=lat_rad*180.0/math.pi
matrix_Ry= np.array([[math.sin(lat_rad), 0,-math.cos(lat_rad)], 
                   [0.0,1.0,0.0],
                   [math.cos(lat_rad),0,math.sin(lat_rad)]]) 
matrix_Rz= np.array([[math.cos(lon_rad),math.sin(lon_rad),0],
                   [-math.sin(lon_rad), math.cos(lon_rad),0],
                   [0.0,0.0,1.0]])
r_sez = np.dot(np.dot(matrix_Ry, matrix_Rz), ecef_vector)
s_km=r_sez[0]
e_km=r_sez[1]
z_km=r_sez[2]
print(s_km)
print(e_km)
print(z_km)