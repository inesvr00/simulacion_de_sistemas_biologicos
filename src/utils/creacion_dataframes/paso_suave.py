def añadir_perdida_energia_paso_suave(df, w, energia):
    df.loc[len(df)] = [w / energia, energia]
    return df