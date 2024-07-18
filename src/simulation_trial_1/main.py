import datetime
import os
import time
from loguru import logger
from matplotlib import pyplot as plt
import pandas as pd
from src.simulation_trial_1.comprobacion_perdida_energia_paso_suave import comprobacion_energia_paso_suave, plot_histogram
from src.simulation_trial_1.kerma_dosis import calcular_rango_efectivo, plot_histogram_kerma_dosis, plot_kerma_and_dosis
from src.simulation_trial_1.photon_simulation import PhotonProperties


    
if __name__ == '__main__':
    start_time = time.time()
    num_photons = 1000000 # Para el análisis de 1,000,000 de fotones introducir este valor
    E_0 = 500
    df_energia_suave =  pd.DataFrame(columns=['r', 'energy'])
    df_dosis =  pd.DataFrame(columns=['energy', 'z'])
    df_kerma = pd.DataFrame(columns=['energy', 'z'])
    df_compton_photo = pd.DataFrame(columns=['compton', 'photo'])
    df_electron_generated = pd.DataFrame(columns=['n_electron', 'absorbed_when_generated'])
    df_elastic_inelastic = pd.DataFrame(columns=['elastic', 'inelastic'])
    step = 0
    for i in range(num_photons):
        step += 1
        photon = PhotonProperties(E_0)  # Crear una nueva instancia para cada fotón.
        df_new_energia_suave, df_new_dosis, df_new_kerma, n_compton, n_photo, n_electron, electron_absorbed, df_new_elastic_inelastic = photon.photon_simulation()
        df_energia_suave = pd.concat([df_energia_suave, df_new_energia_suave], ignore_index=True)
        df_dosis = pd.concat([df_dosis, df_new_dosis], ignore_index=True)
        df_kerma = pd.concat([df_kerma, df_new_kerma], ignore_index=True)
        df_compton_photo.loc[i] = [n_compton, n_photo]
        df_electron_generated.loc[i] = [n_electron, electron_absorbed]
        df_elastic_inelastic = pd.concat([df_elastic_inelastic, df_new_elastic_inelastic], ignore_index=True)
        logger.info(step)

    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    df_energia_suave.to_csv(f'out/csv/energia_suave_{current_time}.csv', index=False)
    df_dosis.to_csv(f'out/csv/dosis_{current_time}.csv', index=False)
    df_kerma.to_csv(f'out/csv/kerma_{current_time}.csv', index=False)
    df_compton_photo.to_csv(f'out/csv/compton_photo{current_time}.csv', index=False)
    df_electron_generated.to_csv(f'out/csv/electron_generated{current_time}.csv', index=False)
    df_elastic_inelastic.to_csv(f'out/csv/elastic_inelastic{current_time}.csv', index=False)
    
    ## Paso correspondiente a la unión de csvs. Debe ser comentado para llevar a cabo una ejecución estándar
    # folder_path = '/home/ines/simulacion_de_sistemas_biologicos/out/csv'
    # csv_files_dosis = [file for file in os.listdir(folder_path) if file.startswith('kerma') and file.endswith('.csv')]
    # dfs = []
    # for file in csv_files_dosis:
    #     file_path = os.path.join(folder_path, file)
    #     df = pd.read_csv(file_path)
    #     dfs.append(df)
    # dosis_total = pd.concat(dfs, ignore_index=True)
    # dosis_total.to_csv('/home/ines/simulacion_de_sistemas_biologicos/out/csv/kerma_total.csv', index=False)
    
    # Obtención de valores de interés a través de los csvs ya concatenados
    # path_compton_photo = '/home/ines/simulacion_de_sistemas_biologicos/out/csv/compton_photo2024-05-16_11-30-24.csv'
    data_compton_photo = pd.read_csv(path_compton_photo)
    suma_compton = data_compton_photo['compton'].sum()
    suma_photo = data_compton_photo['photo'].sum()
    media_compton = data_compton_photo['compton'].mean()
    sem_compton = data_compton_photo['compton'].sem()
    media_photo = data_compton_photo['photo'].mean()
    sem_photo = data_compton_photo['photo'].sem()
    print("Suma de la columna 'compton':", suma_compton)
    print("Suma de la columna 'photo':", suma_photo)
    print("Media por fotón 'compton':", media_compton)
    print("Error estándar de la media para 'compton':", sem_compton)
    print("Media por fotón 'photo':", media_photo)
    print("Error estándar de la media para 'photo':", sem_photo)
    
    # path_electron = '/home/ines/simulacion_de_sistemas_biologicos/out/csv/electron_generated2024-05-16_11-30-24.csv'
    data_electron = pd.read_csv(path_electron)
    suma_electron = data_electron['n_electron'].sum()
    media_electron = data_electron['n_electron'].mean()
    sem_electron = data_electron['n_electron'].sem()
    data_electron['n_e_absorbed_when_generated'] = data_electron['absorbed_when_generated'] * data_electron['n_electron']
    print("Suma de los electrones secundarios generados: ", suma_electron)
    print("Media electrones por fotón:", media_electron)
    print("Error estándar de la media para 'n_electron':", sem_electron)
    print("Porcentaje de los electrones secundarios absorbidos: ", data_electron['n_e_absorbed_when_generated'].sum()/ suma_electron)
    # path_elastic_inelastic = '/home/ines/simulacion_de_sistemas_biologicos/out/csv/elastic_inelastic2024-05-16_11-30-24.csv'
    data_elasctic_inelastic = pd.read_csv(path_elastic_inelastic)
    media_elastic = data_elasctic_inelastic['elastic'].mean()
    sem_elastic = data_elasctic_inelastic['elastic'].sem()
    media_inelastic = data_elasctic_inelastic['inelastic'].mean()
    sem_inelastic = data_elasctic_inelastic['inelastic'].sem()
    print("Media colisiones elásticas por electrón:", media_elastic)
    print("Error estándar de la media colisiones elásticas por electrón:", sem_elastic)
    print("Media colisiones inelásticas por electrón:", media_inelastic)
    print("Error estándar de la media colisiones inelásticas por electrón:", sem_inelastic)
    
    
    
    # df_energia_suave = pd.read_csv("/home/ines/simulacion_de_sistemas_biologicos/out/csv/energia_suave_2024-05-16_11-30-24.csv")
    # df_kerma = pd.read_csv("/home/ines/simulacion_de_sistemas_biologicos/out/csv/kerma_2024-05-16_11-30-24.csv")
    # df_dosis = pd.read_csv("/home/ines/simulacion_de_sistemas_biologicos/out/csv/dosis_2024-05-16_11-30-24.csv")
    # df_kerma_05 = pd.read_csv("/home/ines/simulacion_de_sistemas_biologicos/out/csv/kerma_2024-05-24_00-09-52.csv")
    # df_dosis_05 = pd.read_csv("/home/ines/simulacion_de_sistemas_biologicos/out/csv/dosis_2024-05-24_00-09-52.csv")
    comprobacion_energia_paso_suave(df_energia_suave, num_photons)
    plot_histogram(df_energia_suave, "r")
    plot_histogram_kerma_dosis(df_kerma, "kerma_0_10", num_photons, E_0, 0, 10)
    plot_histogram_kerma_dosis(df_dosis, "dosis_0_10", num_photons, E_0, 0, 10)
    plot_kerma_and_dosis(df_dosis, df_kerma, num_photons, E_0, 0, 10)
    plot_kerma_and_dosis(df_dosis, df_kerma, num_photons, E_0, 0, 100)
    plot_kerma_and_dosis(df_dosis_05, df_kerma_05, num_photons, E_0, 0, 0.5)
    plot_histogram_kerma_dosis(df_kerma, "kerma_0_05", num_photons, E_0, 0, 0.5)
    plot_histogram_kerma_dosis(df_dosis, "dosis_0_05", num_photons, E_0, 0, 0.5)
    plot_histogram_kerma_dosis(df_kerma, "kerma", num_photons, E_0, 0, 500)
    plot_histogram_kerma_dosis(df_dosis, "dosis", num_photons, E_0, 0, 500)
    
    rango_efectivo, z_cruzamiento1, z_cruzamiento2 = calcular_rango_efectivo(df_kerma, df_dosis)

    print(f"Rango efectivo de los electrones (Re): {rango_efectivo}")
    print(f"Primera intersección (z1): {z_cruzamiento1}")
    print(f"Segunda intersección (z2): {z_cruzamiento2}")
    end_time = time.time() 
     
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time}")
