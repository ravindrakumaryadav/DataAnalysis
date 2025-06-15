import numpy as np

def geographic_to_solar_magnetic(lon, lat):
  """Converts from geographic coordinates to solar magnetic coordinates.

  Args:
    lon: The geographic longitude in degrees.
    lat: The geographic latitude in degrees.

  Returns:
    A tuple containing the solar magnetic coordinates (X, Y, Z) in km.
  """

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


