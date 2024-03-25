import math
import random
import numpy as np
from src.simulation.electron_simulation import ElectronProperties

from src.simulation.medium import AbsorbentMedium
from src.simulation.space_geometry import InteractionVolume


class PhotonProperties:
    def __init__(self, E0=500, position_0 = (0, 0, 0), direction_0 = (0, 0, 1), sigma_E0=1.3636e-28, medium = AbsorbentMedium(), geometry = InteractionVolume()):
        """
        Initialize the PhotonProperties with the properties of a photon.

        :param E0: Initial energy of the photon in keV.
        :param omega_E0: Photoelectric cross-section in cm^2.
        """
        self.energy_0 = E0
        self.sigma_E0 = sigma_E0
        self.position_0 = position_0
        self.direction_0 = direction_0
        self.energy = E0
        self.position = position_0
        self.direction = direction_0
        self.step = 0
        self.medium = medium
        self.geometry = geometry
    
    def __str__(self):
        return f"Photon(Energy={self.energy:.2f} keV, Position={self.position}, Direction={self.direction}, Step={self.step})"
        
    def attenuation_coefficient_calculation(self):
        N = self.medium.number_density()
        mu_T = N * (self.effective_section_photoelectric() + self.effective_section_compton())

        return mu_T
        
    def effective_section_photoelectric(self):
        sigma_foto = self.sigma_E0 * (self.energy / self.energy_0)**(-3)

        return sigma_foto
    
    def effective_section_compton(self):
        """
        Calculate the effective Compton atomic cross-section.
        
        :param energy: Energy of the photon in keV.
        :return: Effective Compton atomic cross-section in cm^2.
        """
        mc2 = 510.999  # rest mass energy of electron in keV
        kappa = self.energy / mc2
        r_e = 2.8*10**(-13) # cm
        # The formula for the effective cross-section from the screenshot
        sigma_compton = np.pi * (r_e ** 2) * (2 * kappa * (2 + kappa*(1 + kappa) * (8 + kappa)) + (1 + 2 * kappa)**2 * (kappa * (kappa -2) -2) * math.log(1 + 2 * kappa))/ (kappa**3 * (1 + 2 * kappa)**2)

        return sigma_compton
        
    def free_way_until_next_interaction(self, mu_T):
        """
        Calculate the distance to the next interaction for the photon.
        
        :param mu_T: Total linear attenuation coefficient at the photon's energy.
        :return: Distance to the next interaction.
        """
        U = random.uniform(0, 1)  # Generate a random number U between 0 and 1
        s = -math.log(U) / mu_T

        return s, U
    
    def simulate_scattering_event(self, mu_T,):
        """
        Simulate a scattering event and determine the type of interaction.

        :param energy: Energy of the photon in keV.
        :return: Type of interaction ('photoelectric' or 'compton').
        """    
        N = self.medium.number_density()
        p_photo = N * self.effective_section_photoelectric() / mu_T
        p_co = N * self.effective_section_compton() / mu_T
        
        return p_photo, p_co
        
    def move_photon(self, s):
        """
        Move the photon a distance s in its current direction.
        
        :param s: Distance to move the photon.
        :return: A boolean indicating if the photon has left the interaction volume.
        """
        self.position = tuple(np.array(self.position) + s * np.array(self.direction))

        return self.position[2] > 0  # Returns False if the photon has left the interaction volume
    
    def simulate_photoelectron_emission(self, photoelectron_energy):
        mc2 = 510.999  # rest mass energy of electron in keV
        beta = math.sqrt(photoelectron_energy * (photoelectron_energy + 2 * mc2)) / (photoelectron_energy + mc2)
        A = 1 / beta - 1
        gamma = 1 + photoelectron_energy / mc2

        def g(v):
            return (2 - v) * (1 / (A + v) + 1 / 2 + beta * gamma * (v - 1) * (gamma - 2))

        # Sampling v
        
        while True:
            U1 = random.uniform(0, 1)
            argument = ((A + 2)**2 - 4 * U1 * (A + 2) * U1)

            if argument >= 0:  # Verifica si el argumento de la raíz cuadrada es no negativo
                v = (2 * A / argument)**0.5  # Calcula v solo si es seguro que no resultará en un número imaginario
                break 
        # Sampling U2 and checking g(v)
        U2 = random.uniform(0, 1)
        while U2 >= g(v):
            U1 = random.uniform(0, 1)

            v = (2 * A / ((A + 2)**2 - 4 * U1 * (A + 2) * U1))**0.5

            U2 = random.uniform(0, 1)

        
        theta_e = math.acos(1 - v)
        
        # Sampling azimuthal emission angle phi_e
        phi_e = 2 * math.pi * random.uniform(0, 1)
        
        # New direction vector for the photoelectron
        geometry = InteractionVolume()
        direction = geometry.rotate_vector(self.direction_0, theta_e, phi_e)
        return direction
        
    def simulate_compton_scattering(self):
        phi_p = 2 * math.pi * random.uniform(0, 1)
        phi_e = phi_p + math.pi
        mc2 = 510.999  # rest mass energy of electron in keV
        kappa = self.energy / mc2

        def calculate_tau_min(kappa):
            return 1 / (1 + 2 * kappa)

        def calculate_a1_a2(kappa):
            a1 = np.log(1 + 2 * kappa)
            a2 = 2 * kappa * (1 + kappa) / (1 + 2 * kappa)**2
            return a1, a2

        def select_function(a1, a2):
            U1 = np.random.uniform(0, 1)
            if U1 <= a1 / (a1 + a2):
                return 1  # Select g1
            else:
                return 2  # Select g2

        def sample_tau(tau_min, i):
            U2 = np.random.uniform(0, 1)
            if i == 1:
                return tau_min + U2 * (1 - tau_min)
            else:
                return np.sqrt(tau_min**2 + U2 * (1 - tau_min**2))

        def T_function(tau, kappa):
            return 1 - (1 - tau) * (2 * kappa + 1) / (kappa**2 * (1 + tau**2))

        def calculate_theta(kappa):
            tau_min = calculate_tau_min(kappa)
            a1, a2 = calculate_a1_a2(kappa)

            while True:
                i = select_function(a1, a2)
                tau = sample_tau(tau_min, i)
                T_tau = T_function(tau, kappa)

                U3 = np.random.uniform(0, 1)
                if U3 <= T_tau:
                    cos_theta = 1 - (1 - tau) / (kappa * tau)
                    theta = np.arccos(cos_theta)
                    return theta, cos_theta
        theta, cos_theta = calculate_theta(kappa)

        def calculate_electron_energy(kappa, cos_theta):
            h = 4.135667696e-18  # Planck constant in keV·s
            c = 29979245800  # Speed of light in cm/s
            return (h* c * kappa * (1 - cos_theta)) / (1 + kappa * (1 - cos_theta))

        def calculate_photon_energy(kappa, cos_theta):
            return self.energy / (1 + kappa * (1 - cos_theta))

        def calculate_electron_angle(photon_energy, cos_theta):
            numerator = self.energy - photon_energy * cos_theta
            denominator = np.sqrt(self.energy**2 + photon_energy**2 - 2 * self.energy * photon_energy * cos_theta)
            return np.arccos(numerator / denominator)
        
        photon_energy = calculate_photon_energy(kappa, cos_theta )        
        electron_energy = calculate_electron_energy(kappa, cos_theta)
        theta_electron = calculate_electron_angle(photon_energy, cos_theta)
        print(photon_energy, electron_energy, theta_electron)
        return photon_energy, electron_energy, theta_electron, theta, phi_p, phi_e
        
    def photon_simulation(self):
        continue_simulation = True
        while continue_simulation:

            mu_T = self.attenuation_coefficient_calculation()
            s, U = self.free_way_until_next_interaction(mu_T)
            self.position = self.move_photon(s)
            if self.position:
                p_photo, p_co = self.simulate_scattering_event(mu_T)
                if U <= p_photo: # ES <=
                    ionized_shell = max((shell_energy for shell_energy in [self.medium.U_K, self.medium.U_L1, self.medium.U_L2] if self.energy >= shell_energy), default=0)
                    if ionized_shell > 0:
                        photoelectron_energy = self.energy - ionized_shell
                    else:
                        print("This scenario is very unlikely.")
                        photoelectron_energy = self.energy
                    if photoelectron_energy > self.medium.E_abs_el:
                        photoelectron_direction = self.simulate_photoelectron_emission(photoelectron_energy)
                        electron = ElectronProperties(electron_energy, self.position, photoelectron_direction)
                        continue_simulation = False
                        return electron # Absorbido
                    else:
                        electron = ElectronProperties(self.energy, self.position, self.direction)
                        print(f"The photon has been absorbed at the emission point.")
                        continue_simulation = False
                else:
                    photon_energy, electron_energy, theta_e, theta_p, phi_p, phi_e = self.simulate_compton_scattering()
                    if electron_energy > self.medium.E_abs_el:
                        electron_direction = self.geometry.rotate_vector(self.direction_0, theta_e, phi_e)
                        # Cuidado con las energías, ahora parece que hay que almacenarlas en otros estados iniciales... Tal vez un nuevo objeto de la clase electron y otro de la clase fotón?
                        electron = ElectronProperties(electron_energy, self.position, electron_direction)

                    else:
                        electron = ElectronProperties(electron_energy, self.position, self.direction)

                    if photon_energy > self.medium.E_abs_fo:
                        self.energy = photon_energy
                        self.direction = electron_direction = self.geometry.rotate_vector(self.direction, theta_p, phi_p)       
                        print("Hola")
                    else:
                        self.energy = photon_energy
                        continue_simulation = False
            else:
                print(f'Photon exited the permitted space at position {self.position}...')
        print(electron, self)
        return electron, self    
    
            # return {'new_position': self.position,
            #         'photon_within_volumne': photon_within_volume,
            #         'step': self.step
            # }
            
            
