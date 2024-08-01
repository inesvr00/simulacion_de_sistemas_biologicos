from src.simulation_trial_1.comprobacion_perdida_energia_paso_suave import añadir_perdida_energia_paso_suave
from src.utils.comunes.calcular_angulos import angulo_azimutal
from src.utils.comunes.camino_libre import libre_hasta_siguiente_interaccion
from src.utils.comunes.funciones_espaciales import rotar_vector_direccion
from src.utils.comunes.movimiento_direccion import mover
from src.utils.creacion_dataframes.kerma_dosis import añadir_dosis
from src.utils.electron.calcular_angulos import calcular_theta_elastico
from src.utils.electron.coeficientes_atenuacion import calculo_coeficiente_atenuacion, seccion_eficaz_el, seccion_eficaz_in
from src.utils.electron.colision_inelastica import calcular_energia_perdida_inelastica
from src.utils.electron.paso_suave import calcular_perdida_energia, deflexion_angular
from src.propiedades.electron import propiedades_electron
from src.utils.electron.probabilidades_colision import calcular_p_el, calcular_p_in
import pandas as pd


def historia_electron(posicion, direccion, energia):
    paso = 0
    n_elastico = 0
    n_inelastico = 0
    
    df_energia_suave = pd.DataFrame(columns=['r', 'energia'])
    df_dosis =  pd.DataFrame(columns=['energia', 'z'])
        
    continuar_simulacion = True

    while continuar_simulacion:
        paso += 1
        sigma_el_h = seccion_eficaz_el(energia)
        # print(f"sigma_el_h: {sigma_el_h}")
        sigma_in_h = seccion_eficaz_in(energia)
        # print(f"sigma_in_h: {sigma_in_h}")
        mu_t = calculo_coeficiente_atenuacion(sigma_el_h, sigma_in_h)
        # print(f"mu_t: {mu_t}")
        s, U = libre_hasta_siguiente_interaccion(mu_t)
        # print(f"s: {s}, U: {U}")
        posicion, tau = mover(posicion, direccion, s, 'electron')
        # print(f"posición: {posicion}, tau: {tau}")
        if posicion[2] < 0:
            continuar_simulacion = False
            # print(f"El electrón ha salido del volumen de estudio en {posicion}")
        else:
            chi, phi_may = deflexion_angular(energia, s)
            # print(f"chi: {chi}, phi_may: {phi_may}")
            w = calcular_perdida_energia(energia, s)
            # print(f"w: {w}")

            df_energia_suave = añadir_perdida_energia_paso_suave(df_energia_suave, w, energia)

            direccion = rotar_vector_direccion(direccion, chi, phi_may)
            # print(f"direccion: {direccion}")
            
            df_dosis = añadir_dosis(df_dosis, w, posicion[2])

            # Se deposita w localmente en la posicion
            energia = energia - w
            # print(f"energia: {energia}")
            
            if energia < propiedades_electron['energia_abs']:
                continuar_simulacion = False
                df_dosis = añadir_dosis(df_dosis, energia, posicion[2])
                # La energía se deposita localmente en la posición
            else:
                posicion, valor = mover(posicion, direccion, s-tau, None)
                # print(f"posicion: {posicion}")

                if posicion[2] < 0:
                    continuar_simulacion = False
                    # La partícula ha abandonado el volumen de estudio
                else:
                    p_el = calcular_p_el(sigma_el_h, mu_t)
                    p_in = calcular_p_in(sigma_in_h, mu_t)
                    # print(f"p_el: {p_el}, p_in:{p_in}")
                    
                    if U <= p_el:
                        # print("Elástico")
                        n_elastico += 1
                        theta = calcular_theta_elastico(energia)
                        # print(f"theta: {theta}")
                        phi_min = angulo_azimutal()
                        direccion = rotar_vector_direccion(direccion, theta, phi_min)
                        # print(f"direccion: {direccion}")
                    else:
                        # print(f"Inelástico")
                        n_inelastico += 1
                        W, theta = calcular_energia_perdida_inelastica(energia)
                        # print(f"W: {W}, theta: {theta}")
                        phi_min = angulo_azimutal()
                        direccion = rotar_vector_direccion(direccion, theta, phi_min)
                        # print(f"direccion: {direccion}")
                        if W < energia:
                            energia = energia - W
                            df_dosis = añadir_dosis(df_dosis, W, posicion[2])
                        else:
                            df_dosis = añadir_dosis(df_dosis, energia, posicion[2])
                            energia = 0
                            
                        
                        # W se deposita en el medio en posición
                        if energia < propiedades_electron['energia_abs']:
                            continuar_simulacion = False
                            df_dosis = añadir_dosis(df_dosis, energia, posicion[2])
                            # La energía se deposita en el medio en la posición
    # print(f"Número de choques elásticos: {n_elastico}")
    # print(f"Número de choques inelásticos: {n_inelastico}")
    return df_energia_suave, df_dosis, n_elastico, n_inelastico              
                           
                