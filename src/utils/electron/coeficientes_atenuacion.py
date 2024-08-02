from src.propiedades.medio import propiedades_medio
from src.propiedades.electron import propiedades_electron
import math

def calculo_coeficiente_atenuacion(sigma_el_h, sigma_in_h ):
    """
    Calcula el coeficiente de atenuación total.
    :param sigma_el_h: Sección eficaz efectiva para la dispersión elástica.
    :param sigma_in_h: Sección eficaz efectiva para la dispersión inelástica.
    :return: Coeficiente de atenuación total.
    """
    N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23
    mu_T = N * (sigma_el_h + sigma_in_h)
    return mu_T
    
def seccion_eficaz_el(energia):
    """
    Calcula la sección eficaz efectiva para la dispersión elástica.
    :param energia: Energía del fotón en keV.
    :return: Sección eficaz efectiva para la dispersión elástica.
    
    """
    sigma_el = (1.290e-14) / (1000 * energia)**0.730
    A_0 = 2.5 * (math.log(1000 * energia))**4 / (1000 * energia)**1.434
    mu_c = (1 - math.cos(propiedades_electron['theta_c'])) / 2
    
    sigma_el_h = sigma_el * A_0 * ((1 - mu_c) / (mu_c + A_0))
    return sigma_el_h

def seccion_eficaz_in(energia):
    """
    Calcula la sección eficaz efectiva para la dispersión inelástica.
    :param energia: Energía del fotón en keV.
    :return: Sección eficaz efectiva para la dispersión inelástica.
    """
    mc2 = 510.999
    re = 2.81794e-13
    
    a = (energia / (energia + mc2))**2 
    beta = math.sqrt(energia * (energia + 2 * mc2) / (energia + mc2)**2)
    
    def integrand(W):
        term1 = -(1 / W)
        term2 = a * W / energia**2
        term3 = 1 / (energia - W)
        term4 = ((a - 1) / energia) * math.log(W / (energia - W))
        return (term1 + term2 + term3 + term4)
    
    integrand_value =  integrand(energia / 2) - integrand(propiedades_electron['W_ci'])
    sigma_in_h = propiedades_medio['numero_atomico'] * 2 * math.pi * re**2 * mc2 / beta**2 * integrand_value
    return sigma_in_h     