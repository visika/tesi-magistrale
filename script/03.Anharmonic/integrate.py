"""
Questo script integra con la regola dei trapezoidi i valori della media
dell'energia di Helmholtz sui valori di lambda.

Si aspetta un file di input, averages, con due colonne:
- lambda
- valore della media
"""

import numpy as np

data = np.loadtxt("averages")

lambdas = data[:, 0]
averages = data[:, 1]
# deviations = data[:,2]

dx = 1 / len(lambdas)

integral = dx * (0.5 * averages[0] + averages[1] + averages[2] + 0.5 * averages[3])

print(integral)

# error = np.sqrt( 0.5**2 * ( 0.25 * deviations[0]**2 + deviations[1]**2 + 0.25 * deviations[2]**2 ))
