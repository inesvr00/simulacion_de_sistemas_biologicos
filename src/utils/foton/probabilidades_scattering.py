from src.properties.medio import propiedades_medio


def simular_evento_de_dispersion(seccion_eficaz_fotoelectrica, seccion_eficaz_compton, mu_T):
    """
    Simula un evento de dispersión y determina el tipo de interacción.
    :param seccion_eficaz_fotoelectrica: Sección eficaz fotoeléctrica.
    :param seccion_eficaz_compton: Sección eficaz Compton.
    :param mu_T: Coeficiente de atenuación total a la energía del fotón.
    :return: Probabilidades de cada tipo de interacción ('fotoelectrica' o 'compton').
    """
    N = (propiedades_medio['densidad'] / propiedades_medio['peso_atomico']) * 6.022e23
    p_photo = N * seccion_eficaz_fotoelectrica / mu_T
    p_co = N * seccion_eficaz_compton / mu_T

    return p_photo, p_co
