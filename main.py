from src.simulation.photon_simulation import PhotonProperties


if __name__ == '__main__':
    num_photons = 1
    for i in range(num_photons):
        photon = PhotonProperties()  # Crear una nueva instancia para cada fotón.
        electron, photon = photon.photon_simulation()