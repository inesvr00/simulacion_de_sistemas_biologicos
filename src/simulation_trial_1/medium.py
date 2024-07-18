class AbsorbentMedium:
    def __init__(self, atomic_number=6, density=1.7, atomic_weight=12.0107, 
                 U_K=0.288, U_L1=0.01659, U_L2=0.01126, I=0.078, 
                 E_abs_el=10, E_abs_fo=1):
        """
        Initialize the AbsorbentMedium with the properties of a material.

        :param atomic_number: Atomic number of the medium (default Carbon)
        :param density: Density of the medium in g/cm^3
        :param atomic_weight: Atomic weight of the medium in g/mol
        :param U_K: Ionization energy of the K shell in keV
        :param U_L1: Ionization energy of the L1 shell in keV
        :param U_L2: Ionization energy of the L2 shell in keV
        :param I: Average excitation energy in keV
        :param E_abs_el: Electron absorption energy in keV
        :param E_abs_fo: Photon absorption energy in keV
        """
        self.atomic_number = atomic_number
        self.density = density
        self.atomic_weight = atomic_weight
        self.U_K = U_K  
        self.U_L1 = U_L1  
        self.U_L2 = U_L2  
        self.I = I  
        self.E_abs_el = E_abs_el 
        self.E_abs_fo = E_abs_fo  
        
    def number_density(self):
        """
        Calculate the number density of atoms in the medium.

        :return: Number density of atoms in atoms/cm^3.
        """
        avogadro_number = 6.022e23  # atoms/mol
        number_density = (self.density / self.atomic_weight) * avogadro_number
        return number_density