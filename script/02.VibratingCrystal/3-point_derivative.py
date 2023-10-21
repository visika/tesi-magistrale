"""
Calculates the 3-point-wise derivative, as described in the article by Wilson.
input_name: file with temperature in column 0 equilibrium volume in column 1
"""

input_name = "temp+equilibrium_volume"

import numpy as np
from numpy.polynomial import Polynomial

data = np.loadtxt(input_name)
# print(data)
# print(len(data))

result = []

"""Prendiamo i punti a tre alla volta, ne facciamo un fit parabolico"""
for index in range(len(data) - 2):
    # print(index)
    # print(data[index], data[index+1], data[index+2])
    points = np.array([data[index], data[index + 1], data[index + 2]])
    # print("Points:", points)
    # print(points[:,0])
    coefficients = np.polyfit(points[:, 0], points[:, 1], deg=2)
    # print("Coefficients:", coefficients)
    """Bisogna invertire i coefficienti, perché np.polyfit e Polynomial usano
    due ordinamenti opposti."""
    reversed = np.flip(coefficients)
    fit = Polynomial(reversed)
    derivative = fit.deriv()
    # print("Derivative:", derivative)
    for point in points:
        result.append([point[0], derivative(point[0])])

result = np.array(result)
out = np.column_stack([result[:, 0], result[:, 1] * 1e6])

"""Per i punti in corrispondenza dei quali sono state calcolate più
pendenze, tirare fuori la media"""
u, indices = np.unique(out[:, 0], return_inverse=True)

means = []

for index, unique in enumerate(u):
    points = out[np.where(indices == index)]
    mean = np.mean(points, axis=0)
    # print(mean)
    means.append(mean)

means = np.array(means)
# print(means)

np.savetxt(input_name + ".derivative", means, fmt="%.2f")
