def comprobacion_energia_paso_suave(df_energia_suave, N):
    df_energia_suave_filtrado = df_energia_suave.loc[(df_energia_suave['r'] >= 0) & (df_energia_suave['r'] <= 1)]
    print(len(df_energia_suave_filtrado) / len(df_energia_suave))
    print(f"Porcentaje de energía perdida en el paso suave: {df_energia_suave['energia'].sum() / (500 * N)}")