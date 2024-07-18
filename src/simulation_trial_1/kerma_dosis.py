from matplotlib import pyplot as plt
import numpy as np


def add_kerma(df, energy, z):
    df.loc[len(df)] = [energy, z]
    return df

def add_dosis(df, energy, z):
    df.loc[len(df)] = [energy, z]
    return df

def plot_histogram_kerma_dosis(df, name, N=None, E_0=None, min=None, max=None):
    if N:
        df = df[(df['z'] >= min) & (df['z'] <= max)]
        bin_width = (df['z'].max() - df['z'].min()) / 1000
    
    plt.figure(figsize=(10, 6))
    if N:
        plt.hist(df['z'], bins=1000, weights=df['energy'] / (N * E_0 * bin_width), color='skyblue', edgecolor='black')
    else:
        plt.hist(df['z'], bins=1000, weights=df['energy'], color='skyblue', edgecolor='black')
    plt.xlabel('z')
    plt.ylabel('Sum of energy')
    plt.title(f'Histogram of {name}')
    plt.grid(True)
    plt.savefig(f'/home/ines/simulacion_de_sistemas_biologicos/out/{name}.png')
    
def plot_kerma_and_dosis(df_dosis, df_kerma, N, E_0, min, max):
    # Definir el número de bins y calcular el histograma
    num_bins = 1000
    bins = np.linspace(min, max, num_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # Calcular histogramas y normalizar
    dosis_hist, _ = np.histogram(df_dosis['z'], bins=bins, weights=df_dosis['energy'] * num_bins / (N * E_0))
    kerma_hist, _ = np.histogram(df_kerma['z'], bins=bins, weights=df_kerma['energy'] * num_bins / (N * E_0))

    # Realizar ajustes sobre los histogramas
    kerma_coeffs = np.polyfit(bin_centers, kerma_hist, 1)
    kerma_fit = np.poly1d(kerma_coeffs)
    dosis_coeffs = np.polyfit(bin_centers, dosis_hist, 2)
    dosis_fit = np.poly1d(dosis_coeffs)

    # Encontrar el punto de intersección
    intersection_z = np.roots(dosis_fit - kerma_fit)
    real_roots = intersection_z[np.isreal(intersection_z)].real

    # Filtrar para mantener solo las raíces dentro del rango [min, max] y que sean positivas
    valid_roots = [root for root in real_roots if min <= root <= max and root > 0]
    intersection_points = dosis_fit(valid_roots)

    # Visualización
    plt.figure(figsize=(10, 6))
    plt.bar(bin_centers, kerma_hist, width=bins[1] - bins[0], color='gray', edgecolor='gray', label='Kerma Histogram')
    plt.bar(bin_centers, dosis_hist, width=bins[1] - bins[0], color='black', edgecolor='black', label='Dosis Histogram')
    plt.plot(bin_centers, kerma_fit(bin_centers), 'r--', label='Kerma Fit (linear)')
    plt.plot(bin_centers, dosis_fit(bin_centers), 'b--', label='Dosis Fit (quadratic)')

    for root, point in zip(valid_roots, intersection_points):
        plt.plot(root, point, 'go', markersize=10, label=f'Intersection at z={root:.2f} cm')

    plt.xlabel('z/cm')
    plt.ylabel('Sum of energy normalized')
    plt.title('Histogram of Kerma and Dosis with Fits')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'/home/ines/simulacion_de_sistemas_biologicos/out/dosis_kerma_fits.png')
    plt.show()
    
    
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
