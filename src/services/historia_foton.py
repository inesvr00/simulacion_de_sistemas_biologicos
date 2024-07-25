from src.properties.foton import propiedades_iniciales_foton
from src.properties.medio import propiedades_medio
from src.properties.electron import propiedades_electron
from src.services.historia_electron import historia_electron
from src.utils.comunes.calcular_angulos import angulo_azimutal
from src.utils.comunes.camino_libre import libre_hasta_siguiente_interaccion
from src.utils.foton.coeficientes_atenuacion import calcular_coeficiente_atenuacion, seccion_eficaz_compton, seccion_eficaz_fotoelectrico
from src.utils.foton.efecto_fotoelectrico import angulo_polar, capa_seleccionada
from src.utils.comunes.movimiento_direccion import mover
from src.utils.foton.probabilidades_scattering import simular_evento_de_dispersion
from src.utils.foton.scattering_compton import angulos_foton_electron_dispersado, calcular_energia_electron, calcular_energia_foton, calcular_theta_electron
from src.utils.comunes.funciones_espaciales import rotar_vector_direccion


def historia_foton():
    energia_0 = propiedades_iniciales_foton['E0']
    energia = energia_0
    posicion = propiedades_iniciales_foton['posicion_0']
    direccion = propiedades_iniciales_foton['direccion_0']
    sigma_E0 = propiedades_iniciales_foton['sigma_E0']
    
    n_electron = 0
    
    continuar_simulacion = True
    
    while continuar_simulacion:
        # Reiniciar electrón
        energia_electron = None
        theta_electron = None
        phi_electron = None
        direccion_electron = None
        posicion_electron = None
        
        theta_foton = None
        phi_foton = None

        seccion_eficaz_foto = seccion_eficaz_fotoelectrico(sigma_E0, energia, energia_0)
        seccion_eficaz_co = seccion_eficaz_compton(energia)
        print(f"seccion_eficaz_foto: {seccion_eficaz_foto}, seccion_eficaz_compton: {seccion_eficaz_co}")
        mu_T = calcular_coeficiente_atenuacion(seccion_eficaz_foto, seccion_eficaz_co)
        print(f"mu_T: {mu_T}")
        s, U = libre_hasta_siguiente_interaccion(mu_T)
        print(f"Camino libre: {s}cm")
        posicion, value = mover(posicion, direccion, s, 'foton')
        if posicion[2] < 0:
            continuar_simulacion = False
            print(f"El fotón ha salido del volumen de estudio en {posicion}")
        else:
            p_foto, p_co = simular_evento_de_dispersion(seccion_eficaz_foto, seccion_eficaz_co, mu_T)
            print(f"p_foto: {p_foto}, p_co: {p_co}, U: {U}")

            # Absorción fotoeléctrica
            if U <= p_foto:
                print("Efecto fotoeléctrico")
                n_electron += 1
                energia_electron = capa_seleccionada(energia)
                print(f"Energía electrón: {energia_electron}")
                if energia_electron <= propiedades_electron['energia_abs']:
                    # Electrón absorbido, energía electrón depositada en posición
                    print("Electrón absorbido")
                    pass
                else:
                    phi_electron = angulo_azimutal()
                    theta_electron = angulo_polar(energia_electron)
                    direccion_electron = rotar_vector_direccion(direccion, theta_electron, phi_electron)
                    print(f"angulo azimutal: {phi_electron}, ángulo polar: {theta_electron}")
                    electron_data_after = historia_electron(posicion, direccion_electron, energia_electron)
                energia_foton = 0
                continuar_simulacion = False    

            # Scattering Compton
            else:
                print("Scattering Compton")
                n_electron += 1
                theta_foton, cos_theta, phi_electron, phi_foton, kappa = angulos_foton_electron_dispersado(energia)
                energia_electron = calcular_energia_electron(energia, kappa, cos_theta)
                energia_foton = calcular_energia_foton(energia, kappa, cos_theta)
                theta_electron = calcular_theta_electron(energia, energia_foton, cos_theta)
                print(f"angulo azimutal electrón: {phi_electron}, ángulo polar electrón: {theta_electron}")
                print(f"angulo azimutal fotón: {phi_foton}, ángulo polar fotón: {theta_foton}")


                if energia_electron <= propiedades_electron['energia_abs']:
                    # Electrón absorbido, energía electrón depositada en posición
                    print("Electrón absorbido")
                    pass
                else:
                    direccion_electron = rotar_vector_direccion(direccion, theta_electron, phi_electron)
                    electron_data_after = historia_electron(posicion, direccion_electron, energia_electron)
                    
                if energia_foton <= propiedades_iniciales_foton['energia_abs']:
                    print("Fotón absorbido")
                    # Fotón absorbido
                    continuar_simulacion = False
                else:
                    direccion = rotar_vector_direccion(direccion, theta_foton, phi_foton)
            energia = energia_foton
            
        print(f"Número del electrón: {n_electron}")
        print(f"Datos del electrón -> Dirección: {direccion_electron}, Energía: {energia_electron}")
        print(f"Datos del fotón -> Posicón: {posicion}, Dirección: {direccion}, Energía: {energia}")
        print("-------------------------------------")
        if theta_electron:
            
            pass