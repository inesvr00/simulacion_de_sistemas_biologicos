from src.properties.medio import propiedades_medio
import math
import random

def capa_seleccionada(energia):
    if energia >= propiedades_medio['U_K']:
        energia_electron = energia - propiedades_medio['U_K']
    elif energia >= propiedades_medio['U_L1']:
        energia_electron = energia - propiedades_medio['U_L1']
    elif energia >= propiedades_medio['U_L2']:
        energia_electron = energia - propiedades_medio['U_L2']
    else:
        energia_electron = energia
    return energia_electron

def angulo_azimutal():
    return 2 * math.pi * random.uniform(0, 1)

def angulo_polar(energia_electron):
    beta = math.sqrt(energia_electron * (energia_electron + 2 * 510.999)) / (energia_electron + 510.999)
    A = 1 / beta -1
    gamma = 1 + (energia_electron / 510.999)
    valid = False
    
    def g(v):
        return (2 - v) * ((1 / (A + v)) + (1 / 2 + beta * gamma * (gamma - 1)) * (gamma - 2))
    
    while not valid:
        U_1 = random.uniform(0,1)
        v = ((2 * A) / ((A + 2)**2 - 4 * U_1)) * (2*U_1 + (A + 2) * U_1**(1/2))
        U_2 = random.uniform(0,1)

        if U_2 * g(0) > g(v):
            valid = False
        else:
            valid = True

    return math.acos(1 - v)
