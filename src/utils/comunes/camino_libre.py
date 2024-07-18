import random
import math

def libre_hasta_siguiente_interaccion(mu_T):
    """
    Calcula la distancia hasta la próxima interacción para el fotón.
    
    :param mu_T: Coeficiente de atenuación lineal total a la energía del fotón.
    :return: Distancia hasta la próxima interacción.
    """
    U = random.uniform(0, 1)  # Genera un número aleatorio U entre 0 y 1
    s = (-math.log(U)) / mu_T  # Calcula la distancia s usando la transformada de Laplace

    return s, U
