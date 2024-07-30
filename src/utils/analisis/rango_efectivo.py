import pandas as pd
import numpy as np

def calcular_rango_efectivo(kerma_df, dosis_df):
    z_min = 0
    z_max = 200
    
    if z_min is not None:
        kerma_df = kerma_df[kerma_df['z'] >= z_min]
        dosis_df = dosis_df[dosis_df['z'] >= z_min]
    if z_max is not None:
        kerma_df = kerma_df[kerma_df['z'] <= z_max]
        dosis_df = dosis_df[dosis_df['z'] <= z_max]
        
    kerma_df = kerma_df.sort_values(by='z')
    dosis_df = dosis_df.sort_values(by='z')

    z_common = np.union1d(kerma_df['z'].values, dosis_df['z'].values)
    
    kerma_interp = np.interp(z_common, kerma_df['z'], kerma_df['energy'])
    dosis_interp = np.interp(z_common, dosis_df['z'], dosis_df['energy'])

    diferencia = kerma_interp - dosis_interp
    cruzamientos = np.where(np.diff(np.sign(diferencia)))[0]

    if len(cruzamientos) == 0:
        raise ValueError("Las curvas de kerma y dosis no se cruzan en el rango dado.")
    
    z_cruzamiento1 = z_common[cruzamientos[0]]
    z_cruzamiento2 = z_common[cruzamientos[-1]]

    # Calcular el rango efectivo de los electrones
    rango_efectivo = z_cruzamiento2 - z_cruzamiento1

    return rango_efectivo, z_cruzamiento1, z_cruzamiento2