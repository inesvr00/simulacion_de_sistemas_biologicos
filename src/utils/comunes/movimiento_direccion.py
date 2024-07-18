import numpy as np
import random

def mover(posicion, direccion, s, particula):
    """
    Mueve el fotón una distancia s en su dirección actual.
    
    :param posicion: Posición actual del fotón.
    :param direccion: Dirección de movimiento del fotón.
    :param s: Distancia a mover el fotón.
    :return: Nueva posición.
    """
    if particula == 'electron':
        s = random.uniform(0, s)
    
    posicion = tuple(np.array(posicion) + s * np.array(direccion))

    return posicion
