import numpy as np


def añadir_kerma(df, energia, z):
    df.loc[len(df)] = [energia, z]
    return df

def añadir_dosis(df, energia, z):
    # print(f"Columnas de df_dosis: {df.columns.tolist()}")
    # print(f"Longitud de df_dosis: {len(df)}")
    # print(f"Datos a insertar: {[energia, z]}")
    if energia > 0:
        df.loc[len(df)] = [energia, z]
    return df
