import h5py
import matplotlib as mpl
from matplotlib.colors import LogNorm

mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import numpy as np
from astropy.time import Time

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#--------------solar co ordinate
import astropy.units as u
from sunpy.coordinates import frames

from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates import HeliographicCarrington




def PullVar(fname, vID):
    try:
        with h5py.File(fname, 'r') as hf:
            if vID in hf:
                V = hf[vID][()]
                return V
            else:
                print(f'Dataset {vID} not found in the HDF5 file.')
                return None
    except:
        print(f'HDF5 file {fname} not found.')
        return None




def tStep(fname, nStp=0, aID="time", aDef=0.0):
    with h5py.File(fname, 'r') as hf:
        gID = "Step#%d" % (nStp)
        if aID in hf[gID].attrs:
            t = hf[gID].attrs[aID]
        else:
            t = aDef
    return t

#fname = ('testdeltab-kodik0ra-nb4u0imc.h5');
try:
    fname = 'testdeltab-kodik0ra-hzfdzmzo.h5'

# Extract the necessary variables
    dBn = PullVar(fname, '/Step#2/dBn')
    phi = PullVar(fname, "Phi")
    theta = PullVar(fname, "Theta")

    phi_deg = np.degrees(phi[0])
    theta_deg = np.degrees(theta[0])

    mlat = PullVar(fname, "/Step#2/smlat")
    mlon = PullVar(fname, "/Step#2/smlon")

# Correct the shapes of lon, lat, mlat, mlon, and dBn
    lon = phi_deg

    lat = theta_deg

    mlat = mlat[0]
    mlon = mlon[0]
    dBn = dBn[0]

# Calculate colatitude from latitude
    colat_deg = 90 - lat

# Select the northern hemisphere data
    dBn_northern = dBn[:, 180:]
    phi_deg_northern = lon[:, 180:]
    colat_deg_northern = 90 - colat_deg[:, 180:]
    colat_deg_northern = np.flip(colat_deg_northern, axis=0)

# Extract the northern hemisphere data for mlat and mlon
    nhmlat = mlat[180:, :]
    nhmlon = mlon[180:, :] / (180.0 / np.pi)

# Calculate MJD at the specified step using the tStep function
    mjdStp = tStep(fname, nStp=2, aID="MJD", aDef=-np.inf)

# Convert MJD to a datetime object using the Time class from Astropy
    td = Time(mjdStp, format='mjd').datetime

# This will plot the data on a polar plot
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')

# Generate the colormap and normalization
    dBMin = -80.0
    dBMax = 80.0
    cMap = "seismic"
    vdB = kv.genNorm(dBMin, dBMax, doLog=False)

# Plot the pcolormesh with the correct variable names
    cax = ax.pcolormesh(np.radians(phi_deg_northern), np.radians(colat_deg_northern), dBn_northern, cmap=cMap, norm=vdB)

# Set the title and remove tick labels
#ax.set_title(f'dBn Contour Plot (Northern Hemisphere) - Time: {td:%Y-%m-%d %H:%M:%S} UTC')
    ax.set_title(f'dBn Contour Plot (Northern Hemisphere) - Time: {td:%m-%d-%Y %H:%M:%S} UTC')
    ax.set_xticklabels([])
    ax.set_yticklabels([])

# Add a colorbar at the bottom
    cbar = plt.colorbar(cax, orientation='horizontal', pad=0.08)
    cbar.set_label('dBn')

# Label the degrees around the polar plot
    ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'], fontsize=12)

    plt.savefig('test.png', dpi=200)
    plt.show()


    def geographic_to_solar_magnetic(lon, lat):

        # Converts from geographic coordinates to solar magnetic coordinates.

        # Calculate the geocentric solar magnetic (GSM) coordinates.
        gsm_x = np.cos(lon) * np.cos(lat)
        gsm_y = np.sin(lon) * np.cos(lat)
        gsm_z = np.sin(lat)

        # Rotate the GSM coordinates around the Y-axis by 180 degrees.
        r_y = np.array([[-1, 0, 0],
                        [0, 1, 0],
                        [0, 0, -1]])

        sm_x = gsm_x * r_y[0][0] + gsm_y * r_y[1][0] + gsm_z * r_y[2][0]
        sm_y = gsm_x * r_y[0][1] + gsm_y * r_y[1][1] + gsm_z * r_y[2][1]
        sm_z = gsm_x * r_y[0][2] + gsm_y * r_y[1][2] + gsm_z * r_y[2][2]

        return sm_x, sm_y, sm_z


    sm_x, sm_y, sm_z = geographic_to_solar_magnetic(lon, lat)

    print (f"Solar magnetic coordinates: ({str(sm_x)}, {str(sm_y)}, {str(sm_z)}) km")

    sm_mag = np.sqrt(sm_x + sm_y + sm_z ** 2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.scatter(sm_x, sm_y, s=1, c=sm_mag, cmap='jet', norm=LogNorm())
    ax.set_title('Solar Magnetic Field Magnitude')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')

    cbar = plt.colorbar(ax.collections[0], label='Solar Magnetic Field Magnitude (nT)')

    plt.savefig('solar_magnetic_field_magnitude.png')





except Exception as e:
    print(f'An error occurred: {str(e)}')

#------------------------

    def geographic_to_solar_magnetic(lon, lat):

    #Converts from geographic coordinates to solar magnetic coordinates.


    # Calculate the geocentric solar magnetic (GSM) coordinates.
        gsm_x = np.cos(lon) * np.cos(lat)
        gsm_y = np.sin(lon) * np.cos(lat)
        gsm_z = np.sin(lat)

  # Rotate the GSM coordinates around the Y-axis by 180 degrees.
        r_y = np.array([[-1, 0, 0],
                           [0, 1, 0],
                           [0, 0, -1]])

        sm_x = gsm_x * r_y[0][0] + gsm_y * r_y[1][0] + gsm_z * r_y[2][0]
        sm_y = gsm_x * r_y[0][1] + gsm_y * r_y[1][1] + gsm_z * r_y[2][1]
        sm_z = gsm_x * r_y[0][2] + gsm_y * r_y[1][2] + gsm_z * r_y[2][2]

        return sm_x, sm_y, sm_z

#conversion to image of sm
    fname = 'testdeltab-kodik0ra-hzfdzmzo.h5'

    # Extract the necessary variables
    dBn = PullVar(fname, '/Step#2/dBn')
    phi = PullVar(fname, "Phi")
    theta = PullVar(fname, "Theta")

    phi_deg = np.degrees(phi[0])
    theta_deg = np.degrees(theta[0])

    mlat = PullVar(fname, "/Step#2/smlat")
    mlon = PullVar(fname, "/Step#2/smlon")

    # Correct the shapes of lon, lat, mlat, mlon, and dBn
    lon = phi_deg

    lat = theta_deg

    sm_x, sm_y, sm_z = geographic_to_solar_magnetic(lon, lat)

    print (f"Solar magnetic coordinates: ({str(sm_x)}, {str(sm_y)}, {str(sm_z)}) km")

    sm_mag = np.sqrt(sm_x + sm_y + sm_z ** 2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.scatter(sm_x, sm_y, s=1, c=sm_mag, cmap='jet', norm=LogNorm())
    ax.set_title('Solar Magnetic Field Magnitude')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')

    cbar = plt.colorbar(ax.collections[0], label='Solar Magnetic Field Magnitude (nT)')

    plt.savefig('solar_magnetic_field_magnitude.png')



