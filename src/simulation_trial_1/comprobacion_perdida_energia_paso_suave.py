from matplotlib import pyplot as plt
import pandas as pd

def añadir_perdida_energia_paso_suave(df, w, energy):
    df.loc[len(df)] = [w / energy, energy]
    return df
    

def plot_histogram(df, column, bins=100):
    """
    Función para trazar un histograma de una columna de un DataFrame.
    
    Parámetros:
        - df: DataFrame de pandas.
        - column: Nombre de la columna para la cual se quiere hacer el histograma.
        - bins: Número de contenedores para el histograma (opcional, por defecto 10).
    """
    data = df[column]
    
    plt.hist(data, bins=bins)
    
    plt.xlabel(column)
    plt.ylabel('Frecuencia')
    plt.title(f'Histograma de {column}')

    plt.savefig('/home/ines/simulacion_de_sistemas_biologicos/out/histograma_paso_suave.png')
    
def comprobacion_energia_paso_suave(df_energia_suave, N):
    df_energia_suave_filtrado = df_energia_suave.loc[(df_energia_suave['r'] >= 0) & (df_energia_suave['r'] <= 1)]
    print(len(df_energia_suave_filtrado) / len(df_energia_suave))
    print(f"Porcentaje de energía perdida en el paso suave: {df_energia_suave['energy'].sum() / (500 * N)}")
    