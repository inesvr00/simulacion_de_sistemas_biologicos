from matplotlib import pyplot as plt
import numpy as np

def plot_histograma_kerma_dosis(df, nombre, N=None, E_0=None, min=None, max=None):
    if N:
        df = df[(df['z'] >= min) & (df['z'] <= max)]
        ancho_contenedor = (df['z'].max() - df['z'].min()) / 1000
    
    plt.figure(figsize=(10, 6))
    if N:
        plt.hist(df['z'], bins=1000, weights=df['energia'] / (N * E_0 * ancho_contenedor), color='skyblue', edgecolor='black')
    else:
        plt.hist(df['z'], bins=1000, weights=df['energia'], color='skyblue', edgecolor='black')
    plt.xlabel('z')
    plt.ylabel('Suma de la energía')
    plt.title(f'Histograma de {nombre}')
    plt.grid(True)
    plt.savefig(f'out/plots/{nombre}.png')
    
def plot_kerma_y_dosis(df_dosis, df_kerma, N, E_0, min, max):
    # Definir el número de bins y calcular el histograma
    numero_contenedores = 1000
    bins = np.linspace(min, max, numero_contenedores + 1)
    centro_contenedor = 0.5 * (bins[:-1] + bins[1:])

    # Calcular histogramas y normalizar
    dosis_hist, _ = np.histogram(df_dosis['z'], bins=bins, weights=df_dosis['energia'] * numero_contenedores / (N * E_0))
    kerma_hist, _ = np.histogram(df_kerma['z'], bins=bins, weights=df_kerma['energia'] * numero_contenedores / (N * E_0))

    # Realizar ajustes sobre los histogramas
    kerma_coeffs = np.polyfit(centro_contenedor, kerma_hist, 1)
    kerma_fit = np.poly1d(kerma_coeffs)
    dosis_coeffs = np.polyfit(centro_contenedor, dosis_hist, 2)
    dosis_fit = np.poly1d(dosis_coeffs)

    # Encontrar el punto de intersección
    interseccion_z = np.roots(dosis_fit - kerma_fit)
    raices_reales = interseccion_z[np.isreal(interseccion_z)].real

    # Filtrar para mantener solo las raíces dentro del rango [min, max] y que sean positivas
    raices_validas = [raiz for raiz in raices_reales if min <= raiz <= max and raiz > 0]
    puntos_interseccion = dosis_fit(raices_validas)

    # Visualización
    plt.figure(figsize=(10, 6))
    plt.bar(centro_contenedor, kerma_hist, width=bins[1] - bins[0], color='gray', edgecolor='gray', label='Histograma Kerma')
    plt.bar(centro_contenedor, dosis_hist, width=bins[1] - bins[0], color='black', edgecolor='black', label='Histograma Dosis')
    plt.plot(centro_contenedor, kerma_fit(centro_contenedor), 'r--', label='Ajuste Kerma (lineal)')
    plt.plot(centro_contenedor, dosis_fit(centro_contenedor), 'b--', label='Ajuste Dosis (cuadrático)')

    for raiz, point in zip(raices_validas, puntos_interseccion):
        plt.plot(raiz, point, 'go', markersize=10, label=f'Intersección en z={raiz:.2f} cm')

    plt.xlabel('z/cm')
    plt.ylabel('Suma de la energía normalizada')
    plt.title('Histograma del Kerma y la Dosis con ajustes')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'out/plots/dosis_kerma_fits.png')
    plt.show()