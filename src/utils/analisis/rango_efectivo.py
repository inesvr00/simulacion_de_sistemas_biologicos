import pandas as pd
import numpy as np

def calcular_rango_efectivo(kerma_df, dosis_df):
    z_min = 0.15
    z_max = 0.5

    if z_min is not None:
        kerma_df = kerma_df[kerma_df['z'] >= z_min]
        dosis_df = dosis_df[dosis_df['z'] >= z_min]
    if z_max is not None:
        kerma_df = kerma_df[kerma_df['z'] <= z_max]
        dosis_df = dosis_df[dosis_df['z'] <= z_max]

    kerma_df = kerma_df.sort_values(by='z')
    dosis_df = dosis_df.sort_values(by='z')

    z_common = np.union1d(kerma_df['z'].values, dosis_df['z'].values)

    kerma_interp = np.interp(z_common, kerma_df['z'], kerma_df['energia'])
    dosis_interp = np.interp(z_common, dosis_df['z'], dosis_df['energia'])

    diferencia = kerma_interp - dosis_interp

    rango_efectivo = diferencia.mean()

    return rango_efectivo
