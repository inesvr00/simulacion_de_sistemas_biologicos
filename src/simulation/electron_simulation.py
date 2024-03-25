from src.simulation.medium import AbsorbentMedium


class ElectronProperties:
    def __init__(self, E0, position_0, direction_0, medium = AbsorbentMedium()):
        """
        Initialize the PhotonProperties with the properties of a photon.

        :param E0: Initial energy of the photon in keV.
        :param omega_E0: Photoelectric cross-section in cm^2.
        """
        self.energy_0 = E0
        self.position_0 = position_0
        self.direction_0 = direction_0
        self.energy = E0
        self.position = position_0
        self.direction = direction_0
        self.step = 0
        self.medium = medium
        
    def __str__(self):
        return f"Electron(Energy={self.energy:.2f} keV, Position={self.position}, Direction={self.direction})"
