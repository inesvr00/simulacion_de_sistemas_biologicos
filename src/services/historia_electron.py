

from src.utils.comunes.camino_libre import libre_hasta_siguiente_interaccion
from src.utils.comunes.movimiento_direccion import mover
from src.utils.electron.coeficientes_atenuacion import calculo_coeficiente_atenuacion, seccion_eficaz_el, seccion_eficaz_in


def historia_electron(posicion, direccion, energia):
    paso = 0
        
    continuar_simulacion = True
    
    while continuar_simulacion:
        sigma_el_h = seccion_eficaz_el(energia)
        sigma_in_h = seccion_eficaz_in(energia)
        mu_t = calculo_coeficiente_atenuacion(sigma_el_h, sigma_in_h)
        # print(f"mu_t: {mu_t}")
        s, U = libre_hasta_siguiente_interaccion(mu_t)
        posicion, tau = mover(posicion, direccion, s, 'electron')
        if posicion[2] < 0:
            continuar_simulacion = False
            print(f"El electrón ha salido del volumen de estudio en {posicion}")
        else:
            chi, phi = deflexion_angular(s, sigma_el_h)     