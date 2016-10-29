import pandas as pd
'''
Gets data from NOAA's ERDDAP and saves it into .csv file
'''

# Temp, Salt, Oxy
#url = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/nodcWoa09mon1t_LonPM180.csv?temperature_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)],salinity_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)],oxygenSat_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)],oxygenSat_ma[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)],disOxygen_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)]'

# Nitrate
url = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/nodcWoa09mon1n_LonPM180.csv?nitrate_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0.0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)],phosphate_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0.0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)],silicate_an[(0000-01-16):1:(0000-12-16T00:00:00Z)][(0.0):1:(0)][(44.5):1:(44.5)][(-61.5):1:(-61.5)]'
data = pd.read_csv(url)

data.to_csv('Nitrate.csv')