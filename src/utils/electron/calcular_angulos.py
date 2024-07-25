import random
import math


def calcular_theta_elastico(energia):
    """
    Calcular el ángulo de dispersión theta para una colisión elástica.

    :param U: Variable aleatoria.
    :return: Ángulo de dispersión theta.
    """
    U = random.uniform(0, 1)
    A_0 = 2.5 * (math.log10(1000 * energia))**4 / (1000 * energia)**1.434
    return A_0 * U / (1 + A_0 - U)


