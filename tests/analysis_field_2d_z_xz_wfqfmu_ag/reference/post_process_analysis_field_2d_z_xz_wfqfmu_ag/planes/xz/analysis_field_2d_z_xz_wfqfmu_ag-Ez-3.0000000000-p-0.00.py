import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

csv_file = 'analysis_field_2d_z_xz_wfqfmu_ag-Ez-3.0000000000-p-0.00.csv'
output = 'analysis_field_2d_z_xz_wfqfmu_ag-Ez-3.0000000000-p-0.00.svg'
output_png = 'analysis_field_2d_z_xz_wfqfmu_ag-Ez-3.0000000000-p-0.00.png'

data = []
with open(csv_file, 'r') as file:
    for line in file:
        if line.strip():
            data.append([float(line[:25]), float(line[28:53]), float(line[56:81])])
data = np.array(data)

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

point = 500
xi = np.linspace(min(x), max(x), point)
yi = np.linspace(min(y), max(y), point)
X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z, (X, Y), method='linear')

fig = plt.figure(figsize=(10, 10))
plt.pcolormesh(X, Y, Z, cmap='jet', shading='auto')
cbar = plt.colorbar(extend='both')
cbar.formatter.set_useMathText(True)
cbar.update_ticks()
plt.xlabel('x [Å]')
plt.ylabel('z [Å]')
plt.savefig(output_png, dpi=300)
plt.savefig(output, format='svg')
plt.close()
