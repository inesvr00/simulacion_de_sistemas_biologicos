import time
from loguru import logger
import pandas as pd
import datetime

from src.services.historia_foton import historia_foton
from src.propiedades.foton import propiedades_iniciales_foton


def historia_completa_dataframe():
    num_fotones = 100000
    E_0 = propiedades_iniciales_foton['E0']
    df_energia_suave = pd.DataFrame(columns=['r', 'energia'])
    df_dosis = pd.DataFrame(columns=['energia', 'z'])
    df_kerma = pd.DataFrame(columns=['energia', 'z'])
    df_compton_foto = pd.DataFrame(columns=['compton', 'foto'])
    df_electron_generado = pd.DataFrame(columns=['n_electron', 'absorbido_al_generarse'])
    df_elastico_inelastico = pd.DataFrame(columns=['elastico', 'inelastico'])
    paso = 0
    
    for i in range(num_fotones):
        paso += 1
        df_nueva_energia_suave, df_nueva_dosis, df_nueva_kerma, n_compton, n_foto, n_electron, electron_absorbido, df_nuevo_elastico_inelastico = historia_foton()  
        df_energia_suave = pd.concat([df_energia_suave, df_nueva_energia_suave], ignore_index=True)
        df_dosis = pd.concat([df_dosis, df_nueva_dosis], ignore_index=True)
        df_kerma = pd.concat([df_kerma, df_nueva_kerma], ignore_index=True)
        df_compton_foto.loc[i] = [n_compton, n_foto]
        df_electron_generado.loc[i] = [n_electron, electron_absorbido]
        df_elastico_inelastico = pd.concat([df_elastico_inelastico, df_nuevo_elastico_inelastico], ignore_index=True)
        logger.info(paso)
    
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    df_energia_suave.to_csv(f'out/csv/energia_suave_{current_time}.csv', index=False)
    df_dosis.to_csv(f'out/csv/dosis_{current_time}.csv', index=False)
    df_kerma.to_csv(f'out/csv/kerma_{current_time}.csv', index=False)
    df_compton_foto.to_csv(f'out/csv/compton_foto_{current_time}.csv', index=False)
    df_electron_generado.to_csv(f'out/csv/electron_generado_{current_time}.csv', index=False)
    df_elastico_inelastico.to_csv(f'out/csv/elastico_inelastico_{current_time}.csv', index=False)
    