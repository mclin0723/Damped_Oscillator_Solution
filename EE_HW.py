import numpy as np
import matplotlib.pyplot as plt

# Parameters
omega_0 = 1  # Natural frequency
zeta_values = [0.1, 0.5, 1, 2, 3]  # Damping ratios
t = np.linspace(0, 30, 500)  # Time range for omega_0 * t

# Define the solution based on damping ratio
def damped_oscillation(t, omega_0, zeta):
    if zeta < 1:  # Under-damped case
        omega_d = omega_0 * np.sqrt(1 - zeta**2)
        k1 = -1
        k2 = -zeta / np.sqrt(1 - zeta**2)
        return np.exp(-zeta * omega_0 * t) * (k1 * np.cos(omega_d * t) + k2 * np.sin(omega_d * t))
    elif zeta == 1:  # Critically damped case
        k1 = -1
        k2 = 0
        return (k1 + k2 * t) * np.exp(-zeta * omega_0 * t)
    else:  # Over-damped case
        r1 = omega_0 * (-zeta + np.sqrt(zeta**2 - 1))
        r2 = omega_0 * (-zeta - np.sqrt(zeta**2 - 1))
        k1 = (-zeta - np.sqrt(zeta ** 2 - 1)) / (2 * np.sqrt(zeta ** 2 - 1))
        k2 = (zeta + np.sqrt(zeta ** 2 - 1)) / (2 * np.sqrt(zeta ** 2 - 1)) - 1
        # The response is a combination of two exponentials based on initial conditions
        return (k1 * np.exp(r1 * t) + k2 * np.exp(r2 * t))

# Plotting
plt.figure(figsize=(10, 6))

for zeta in zeta_values:
    x_t = damped_oscillation(t, omega_0, zeta) + 1
    plt.plot(t, x_t, label=f'ζ = {zeta}')

plt.title('Damped Oscillator Solution for Different Damping Ratios (ζ)')
plt.xlabel(r'$\omega_0 t$')
plt.ylabel(r'$x(t)$')
plt.legend()
plt.grid(True)
plt.axis([0, 30, 0, 2])
plt.show()