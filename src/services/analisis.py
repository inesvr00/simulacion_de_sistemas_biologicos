import os
import pandas as pd

from src.utils.analisis.paso_suave import comprobacion_energia_paso_suave
from src.utils.analisis.rango_efectivo import calcular_rango_efectivo
from src.utils.plots.kerma_dosis import plot_histograma_kerma_dosis, plot_kerma_y_dosis
from src.utils.plots.paso_suave import plot_histogram


def analisis():
    def obtener_primer_archivo(directorio, prefijo):
        archivos = [f for f in os.listdir(directorio) if f.startswith(prefijo)]
        if archivos:
            return os.path.join(directorio, archivos[0])
        else:
            print(f"No se encontró ningún archivo que comience con '{prefijo}' en el directorio especificado.")
            return None

    directorio_base = "out/csv"
    prefijos = ["compton_foto", "elastico_inelastico", "electron_generado", "dosis", "kerma", "energia_suave"]

    dataframes = {}

    for prefijo in prefijos:
        ruta = obtener_primer_archivo(directorio_base, prefijo)
        if ruta:
            dataframes[prefijo] = pd.read_csv(ruta)
            print(f"Cargado archivo para {prefijo} desde: {ruta}")

    data_compton_foto = dataframes.get("compton_foto")
    data_elastico_inelastico = dataframes.get("elastico_inelastico")
    data_electron = dataframes.get("electron_generado")
    df_dosis = dataframes.get("dosis")
    df_kerma = dataframes.get("kerma")
    df_energia_suave = dataframes.get("energia_suave")
    
    print(data_compton_foto.head())
    print(data_compton_foto.dtypes)
    suma_compton = data_compton_foto['compton'].sum()
    suma_foto = data_compton_foto['foto'].sum()
    media_compton = data_compton_foto['compton'].mean()
    sem_compton = data_compton_foto['compton'].sem()
    media_foto = data_compton_foto['foto'].mean()
    sem_foto = data_compton_foto['foto'].sem()
    print("Suma de la columna 'compton':", suma_compton)
    print("Suma de la columna 'foto':", suma_foto)
    print("Media por fotón 'compton':", media_compton)
    print("Error estándar de la media para 'compton':", sem_compton)
    print("Media por fotón 'foto':", media_foto)
    print("Error estándar de la media para 'foto':", sem_foto)

    suma_electron = data_electron['n_electron'].sum()
    media_electron = data_electron['n_electron'].mean()
    sem_electron = data_electron['n_electron'].sem()
    data_electron['n_e_absorbed_when_generated'] = data_electron['absorbido_al_generarse'] * data_electron['n_electron']
    print("Suma de los electrones secundarios generados: ", suma_electron)
    print("Media electrones por fotón:", media_electron)
    print("Error estándar de la media para 'n_electron':", sem_electron)
    print("Porcentaje de los electrones secundarios absorbidos: ", data_electron['absorbido_al_generarse'].sum()/ suma_electron)

    media_elastico = data_elastico_inelastico['elastico'].mean()
    sem_elastico = data_elastico_inelastico['elastico'].sem()
    media_inelastico = data_elastico_inelastico['inelastico'].mean()
    sem_inelastico = data_elastico_inelastico['inelastico'].sem()
    print("Media colisiones elásticas por electrón:", media_elastico)
    print("Error estándar de la media colisiones elásticas por electrón:", sem_elastico)
    print("Media colisiones inelásticas por electrón:", media_inelastico)
    print("Error estándar de la media colisiones inelásticas por electrón:", sem_inelastico)
    
    n_fotones = 1000000
    E_0 = 500
    comprobacion_energia_paso_suave(df_energia_suave, n_fotones)
    plot_histogram(df_energia_suave, "r")
    plot_histograma_kerma_dosis(df_kerma, "kerma_0_10", n_fotones, E_0, 0, 10)
    plot_histograma_kerma_dosis(df_dosis, "dosis_0_10", n_fotones, E_0, 0, 10)
    plot_kerma_y_dosis(df_dosis, df_kerma, n_fotones, E_0, 0, 10)
    plot_kerma_y_dosis(df_dosis, df_kerma, n_fotones, E_0, 0, 100)
    plot_kerma_y_dosis(df_dosis, df_kerma, n_fotones, E_0, 0, 0.5)
    plot_histograma_kerma_dosis(df_kerma, "kerma_0_05", n_fotones, E_0, 0, 0.5)
    plot_histograma_kerma_dosis(df_dosis, "dosis_0_05", n_fotones, E_0, 0, 0.5)
    plot_histograma_kerma_dosis(df_kerma, "kerma", n_fotones, E_0, 0, 500)
    plot_histograma_kerma_dosis(df_dosis, "dosis", n_fotones, E_0, 0, 500)
    
    rango_efectivo, z_cruzamiento1, z_cruzamiento2 = calcular_rango_efectivo(df_kerma, df_dosis)

    print(f"Rango efectivo de los electrones (Re): {rango_efectivo}")
    print(f"Primera intersección (z1): {z_cruzamiento1}")
    print(f"Segunda intersección (z2): {z_cruzamiento2}")
