import numpy as np
import math
import random


def rotar_vector_direccion(vector_original, theta, phi):
    """
    Rotar un vector por los ángulos dados theta y phi.

    Parámetros:
    vector_original (array): El vector original a rotar.
    theta (float): El ángulo para rotar alrededor del eje x.
    phi (float): El ángulo para rotar alrededor del eje z.

    Retorna:
    numpy.ndarray: El nuevo vector después de la rotación.
    """
    
    # Matrices de rotación
    R_y_minus_theta = np.array([[np.cos(-theta), 0, np.sin(-theta)], 
                                [0, 1, 0], 
                                [-np.sin(-theta), 0, np.cos(-theta)]])
    R_z_minus_phi = np.array([[np.cos(-phi), -np.sin(-phi), 0], 
                             [np.sin(-phi), np.cos(-phi), 0], 
                             [0, 0, 1]])
    
    # Comprobar si el vector es paralelo o antiparalelo al eje z
    if vector_original[2] == 1 :
        u_prima = np.sin(theta) * np.cos(phi)
        v_prima = np.sin(theta) * np.sin(phi)
        w_prima = np.cos(theta)
        
    elif vector_original[2] == -1:
        u_prima = -np.sin(theta) * np.cos(phi)
        v_prima = -np.sin(theta) * np.sin(phi)
        w_prima = -np.cos(theta)
    
    else:
        # Extraer componentes individuales del vector original
        u, v, w = vector_original

        # Calcular los nuevos componentes del vector después de la rotación
        u_prima = u * np.cos(theta) + (u * v * np.cos(phi) - v * np.sin(phi)) * np.sin(theta) / np.sqrt(1 - w**2)
        v_prima = v * np.cos(theta) + (v * w * np.cos(phi) + u * np.sin(phi)) * np.sin(theta) / np.sqrt(1 - w**2)
        w_prima = w * np.cos(theta) - np.sqrt(1 - w**2) * np.sin(theta) * np.cos(phi)

    nueva_direccion = np.array([u_prima, v_prima, w_prima])
    
    # Normalizar el vector para asegurar que tenga módulo unidad
    magnitud = np.linalg.norm(nueva_direccion)
    nueva_direccion_normalizada = nueva_direccion / magnitud

    return nueva_direccion_normalizada



