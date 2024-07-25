from src.properties.medio import propiedades_medio


def calcular_p_el(sigma_el_h, mu_t):
    """
    Calcular la probabilidad de una colisión elástica.

    :param sigma_el_h: Sección eficaz para dispersión elástica.
    :param mu_t: Coeficiente de atenuación total.
    :return: Probabilidad de una colisión elástica.
    """
    
    N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23

    return N * sigma_el_h / mu_t

def calcular_p_in(sigma_in_h, mu_t):
    """
    Calcular la probabilidad de una colisión inelástica.

    :param sigma_in_h: Sección eficaz para dispersión inelástica.
    :param mu_t: Coeficiente de atenuación total.
    :return: Probabilidad de una colisión inelástica.
    """

    N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23
    
    return N * sigma_in_h / mu_t
