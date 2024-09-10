from src.services.analisis import analisis
from src.services.historia_completa_dataframes import historia_completa_dataframe
import warnings


if __name__ == '__main__':
    warnings.simplefilter(action='ignore', category=FutureWarning)

    historia_completa_dataframe()

    analisis()
