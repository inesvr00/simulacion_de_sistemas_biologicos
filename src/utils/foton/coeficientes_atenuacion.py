import numpy as np
import math
from src.propiedades.medio import propiedades_medio

def calcular_coeficiente_atenuacion(seccion_eficaz_fotoelectrica, seccion_eficaz_compton):
    """
    Calcula el coeficiente de atenuación total en base a la sección eficaz
    fotoeléctrica y de compton
    .
    :return: Coeficiente de atenuación lineal total para la energía del fotón.
    """
    N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23
    mu_T = N * (seccion_eficaz_fotoelectrica + seccion_eficaz_compton)
    return mu_T
    
def seccion_eficaz_fotoelectrico(sigma_E0, energia, energia_0):
    """
    Calcular la sección eficaz fotoeléctrica.
    :return: Sección eficaz fotoeléctrica en cm^2.
    """
    sigma_foto = sigma_E0 * (energia / energia_0)**(-3)
    return sigma_foto

def seccion_eficaz_compton(energia):
        """
        Calcular la sección efizac Compton.

        :return: Sección efizac Compton en cm^2.
        """
        kappa = energia / 510.999
        r_e = 2.8*10**(-13) # cm
        
        seccion_eficaz_total_electronica = np.pi * (r_e ** 2) * (2 * kappa * (2 + kappa*(1 + kappa) * (8 + kappa)) + (1 + 2 * kappa)**2 * (kappa * (kappa -2) -2) * math.log(1 + 2 * kappa)) / (kappa**3 * (1 + 2 * kappa)**2)

        return seccion_eficaz_total_electronica * propiedades_medio['numero_atomico']