from src.properties.medio import propiedades_medio
from src.properties.electron import propiedades_electron
import math
import random

def deflexion_angular(energia, s, sigma_el_h):
    """
    Calcular la deflexión angular del electrón después de una colisión.

    :param s: Distancia recorrida por el electrón.
    :param sigma_el_h: Sección eficaz para dispersión elástica.
    :return: Los ángulos de deflexión chi y phi.
    """

    def integrando_deflexion(energia):
        N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23
        A_0 = 2.5 * (math.log10(1000 * energia))**4 / (1000 * energia)**1.434
        mu_c = (1 - math.cos(propiedades_electron['theta_c'])) / 2
        sigma_el = 1.290 * 10**(-14) / (1000 * energia)**0.730
        lambda_el_1 = 1 / (2 * N * sigma_el * A_0 * (1 + A_0) * (math.log((mu_c + A_0) / A_0) - mu_c / (mu_c + A_0)))
        return lambda_el_1

    chi = math.acos(math.exp(-s / integrando_deflexion(energia)))
    phi = 2 * math.pi * random.uniform(0, 1)
    return chi, phi

def calcular_perdida_energia(energia, s):
    """
    Calcular la pérdida de energía del electrón después de recorrer una distancia s.

    :param s: Distancia recorrida por el electrón.
    :return: La energía perdida por el electrón.
    """
    N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23
    Z = propiedades_medio['numero_atomico']
    I = propiedades_medio['I']
    mc2 = 510.999
    re = 2.81794e-13
    
    gamma = (energia + mc2) / mc2
    beta = math.sqrt((gamma**2 - 1) / gamma**2)

    # Corrección de densidad de Fermi 
    delta_f = 0 

    S_s_E = N * Z * 2 * math.pi * re**2 * mc2 / beta**2 * (math.log((energia**2 / I**2) * (gamma + 1) / 2) + 1 - beta**2 - (2 * gamma - 1) / gamma**2 * math.log(2) + 1/8 * ((gamma - 1) / gamma)**2 - delta_f)

    w = S_s_E * s
    
    if w >= energia:
        return energia
    else:
        return w
