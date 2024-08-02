from matplotlib import pyplot as plt

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

    plt.savefig('/home/debian/dev/simulacion_de_sistemas_biologicos/out/plots/histograma_paso_suave.png')