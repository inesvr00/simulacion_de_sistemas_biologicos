import time
from loguru import logger
import pandas as pd
import datetime
import multiprocessing as mp
from src.services.historia_foton import historia_foton
from src.propiedades.foton import propiedades_iniciales_foton
from tqdm import tqdm


def worker(_):
    return historia_foton()


def historia_completa_dataframe():
    num_fotones = 1000000
    num_procesos = mp.cpu_count()
    
    E_0 = propiedades_iniciales_foton['E0']
    df_energia_suave = pd.DataFrame(columns=['r', 'energia'])
    df_dosis = pd.DataFrame(columns=['energia', 'z'])
    df_kerma = pd.DataFrame(columns=['energia', 'z'])
    df_compton_foto = pd.DataFrame(columns=['compton', 'foto'])
    df_electron_generado = pd.DataFrame(columns=['n_electron', 'absorbido_al_generarse'])
    df_elastico_inelastico = pd.DataFrame(columns=['elastico', 'inelastico'])
    paso = 0

    resultados_energia_suave = []
    resultados_dosis = []
    resultados_kerma = []
    resultados_compton_foto = []
    resultados_electron_generado = []
    resultados_elastico_inelastico = []
    
    
    with mp.Pool(num_procesos) as pool:
        # pool.imap es un iterador que permite la actualización de la barra de progreso
        for result in tqdm(pool.imap(worker, [None] * num_fotones), total=num_fotones):
            df_nueva_energia_suave, df_nueva_dosis, df_nueva_kerma, n_compton, n_foto, n_electron, n_electron_absorbido, df_nuevo_elastico_inelastico = result
            
            resultados_energia_suave.append(df_nueva_energia_suave)
            resultados_dosis.append(df_nueva_dosis)
            resultados_kerma.append(df_nueva_kerma)
            resultados_compton_foto.append((n_compton, n_foto))
            resultados_electron_generado.append((n_electron, n_electron_absorbido))
            resultados_elastico_inelastico.append(df_nuevo_elastico_inelastico)
        
    df_energia_suave = pd.concat(resultados_energia_suave, ignore_index=True)
    df_dosis = pd.concat(resultados_dosis, ignore_index=True)
    df_kerma = pd.concat(resultados_kerma, ignore_index=True)
    df_compton_foto = pd.DataFrame(resultados_compton_foto, columns=['compton', 'foto'])
    df_electron_generado = pd.DataFrame(resultados_electron_generado, columns=['n_electron', 'absorbido_al_generarse'])
    df_elastico_inelastico = pd.concat(resultados_elastico_inelastico, ignore_index=True)
    
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    df_energia_suave.to_csv(f'out/csv/energia_suave_{current_time}.csv', index=False)
    df_dosis.to_csv(f'out/csv/dosis_{current_time}.csv', index=False)
    df_kerma.to_csv(f'out/csv/kerma_{current_time}.csv', index=False)
    df_compton_foto.to_csv(f'out/csv/compton_foto_{current_time}.csv', index=False)
    df_electron_generado.to_csv(f'out/csv/electron_generado_{current_time}.csv', index=False)
    df_elastico_inelastico.to_csv(f'out/csv/elastico_inelastico_{current_time}.csv', index=False)
    