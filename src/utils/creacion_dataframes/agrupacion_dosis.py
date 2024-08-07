def agrupar_y_sumar_energia(df):
    df['z_redondeado'] = df['z'].round(3)  # Redondea la columna 'z' a 3 decimales
    df['grupo'] = (df['z_redondeado'].diff().abs() > 0.001).cumsum()  # Crea un grupo para valores consecutivos
    resultado = df.groupby('grupo').agg({'energia': 'sum', 'z_redondeado': 'first'}).reset_index(drop=True)  # Agrupa y suma
    resultado.rename(columns={'z_redondeado': 'z'}, inplace=True)  # Renombra la columna 'z_redondeado' a 'z'
    return resultado[['z', 'energia']]
