import random
import math
import numpy as np
from src.properties.electron import propiedades_electron

def calcular_energia_perdida_inelastica(energia):
    """
    Calcular la energía perdida durante una colisión inelástica.

    :return: Energía perdida y ángulo de dispersión theta.
    """
    mc2 = 510.999
    k_c = propiedades_electron['W_ci'] / energia

    a = (energia / (energia + mc2))**2

    continuar_funcion = True
    p_1 = k_c / (1 - 2 * k_c) * k_c**(-2)
    p_2 = 2 / (1 - 2 * k_c)
    
    while continuar_funcion:
        U_1 = random.uniform(0, 1)
        U_2 = random.uniform(0, 1)

        if U_1 < (1 / (1 + 5 * a * k_c / 2)):
            p_k = p_1
            k = k_c / (1 - U_2 * (1 - 2 * k_c))
        else:
            p_k = p_2
            k = k_c + (U_2 * (1 - 2 * k_c) / 2)

        U_3 = random.uniform(0, 1)
        if U_3 * calcular_phi_k(k_c, k, a) <= calcular_f_in_k(k_c, k, a):
            continuar_funcion = False

    W = k * energia
    theta = np.arccos((energia - W) / energia * (energia + 2 * mc2) / (energia - W + 2 * mc2))

    return W, theta


def calcular_phi_k(k_c, k, a):
    """
    Calcular la función phi_k utilizada en los cálculos de colisiones inelásticas.

    :param k_c: Relación de energía de corte.
    :param k: Relación de energía.
    :param a: Parámetro a.
    :return: Valor de phi_k.
    """

    phi_k = (k**-2 + 5 * a) * np.heaviside(k - k_c, 1) * np.heaviside(1/2 - k, 1)
    return phi_k


def calcular_f_in_k(k_c, k, a):
    """
    Calcular la función f_in_k utilizada en los cálculos de colisiones inelásticas.

    :param k_c: Relación de energía de corte.
    :param k: Relación de energía.
    :param a: Parámetro a.
    :return: Valor de f_in_k.
    """

    f_in_k = ((1 / k**2) + (1 / (1 - k)**2) - (1 / (k * (1 - k))) + a * (1 + 1 / (k * (1 - k)))) * np.heaviside(k - k_c, 1) * np.heaviside(1/2 - k, 1)
    return f_in_k
