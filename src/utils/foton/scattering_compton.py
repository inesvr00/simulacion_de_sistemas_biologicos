import numpy as np
import math
import random


def angulos_foton_electron_dispersado(energia):
    phi_p = 2 * math.pi * random.uniform(0, 1)
    phi_e = phi_p + math.pi

    kappa = energia / 510.999

    def calcular_tau_min(kappa):
        return 1 / (1 + 2 * kappa)

    def calcular_a1_a2(kappa):
        a1 = np.log(1 + 2 * kappa)
        a2 = (2 * kappa * (1 + kappa)) / (1 + 2 * kappa)**2
        return a1, a2

    def seleccionar_funcion(a1, a2):
        U1 = random.uniform(0, 1)
        if 0 <= U1 <= a1 / (a1 + a2):
            return 1  # Selecciona g1
        else:
            return 2  # Selecciona g2

    def muestrear_tau(tau_min, i):
        U2 = random.uniform(0, 1)
        if i == 1:
            return tau_min**U2
        else:
            return np.sqrt(tau_min**2 + U2 * (1 - tau_min**2))

    def funcion_T(tau, kappa):
        return 1 - ((1 - tau) * (tau * (2 * kappa + 1) - 1)) / (kappa**2 * tau * (1 + tau**2))

    def calcular_theta(kappa):
        tau_min = calcular_tau_min(kappa)
        a1, a2 = calcular_a1_a2(kappa)

        while True:
            i = seleccionar_funcion(a1, a2)
            tau = muestrear_tau(tau_min, i)
            T_tau = funcion_T(tau, kappa)

            U3 = random.uniform(0, 1)
            if U3 <= T_tau:
                cos_theta = 1 - ((1 - tau) / (kappa * tau))
                theta = math.acos(cos_theta)
                return theta, cos_theta

    theta, cos_theta = calcular_theta(kappa)
    
    return theta, cos_theta, phi_e, phi_p, kappa

def calcular_energia_electron(energia, kappa, cos_theta):
    return ((energia * kappa * (1 - cos_theta)) / (1 + kappa * (1 - cos_theta)))

def calcular_energia_foton(energia, kappa, cos_theta):
    return (energia / (1 + kappa * (1 - cos_theta)))

def calcular_theta_electron(energia, energia_foton, cos_theta):
    numerator = energia - energia_foton * cos_theta
    denominator = np.sqrt(energia**2 + energia_foton**2 - 2 * energia * energia_foton * cos_theta)
    return np.arccos(numerator / denominator)