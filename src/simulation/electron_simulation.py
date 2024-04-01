import math
import random

import numpy as np
import scipy
from src.simulation.medium import AbsorbentMedium
from src.simulation.space_geometry import InteractionVolume



class ElectronProperties:
    def __init__(self, E0, position_0, direction_0, theta_c=20, W_ci=1, medium = AbsorbentMedium(), geometry = InteractionVolume()):
        """
        Initialize the PhotonProperties with the properties of a photon.

        :param theta_c: Cut angle for elastic collisions.
        :param W_ci: Cut energy for inelastic collisions.
        """
        self.energy_0 = E0
        self.position_0 = position_0
        self.direction_0 = direction_0
        self.energy = E0
        self.position = position_0
        self.direction = direction_0
        self.step = 0
        self.medium = medium
        self.theta_c = theta_c
        self.W_ci = W_ci
        self.A_0 = 2.5 * (math.log10(1000 * self.energy))**4 / (1000 * self.energy)**1.434
        self.mu_c = (1 - math.cos(self.theta_c)) / 2
        self.geometry = geometry
        
    def __str__(self):
        return f"Electron(Energy={self.energy:.2f} keV, Position={self.position}, Direction={self.direction})"

    def attenuation_coefficient_calculation(self, sigma_el_h, sigma_in_h ):
        N = self.medium.number_density()
        mu_T = N * (sigma_el_h + sigma_in_h)

        return mu_T
        
    def effective_section_el(self):
        sigma_el = (1.290e-14) / (1000 * self.energy)**0.730
        sigma_el_h = sigma_el * self.A_0 * (1 - self.mu_c) / (self.mu_c + self.A_0)
        print(f"sigma_el_h: {sigma_el_h}")
        return sigma_el_h
    
    def effective_section_in(self):
        mc2 = 510.999
        re = 2.81794e-13
        
        a = (self.energy / (self.energy + mc2))**2 
        beta = np.sqrt(self.energy * (self.energy + 2 * mc2) / (self.energy + mc2)**2)
        
        def integrand(W):
            term1 = -(1 / W)
            term2 = a * W / self.energy**2
            term3 = 1 / (self.energy - W)
            term4 = ((a - 1) / self.energy) * np.log(W / (self.energy - W))
            return (term1 + term2 + term3 + term4) 
        integrand_value =  integrand(self.energy / 2) - integrand(self.W_ci)

        sigma_in_h = 2 * np.pi * re**2 * mc2 / beta**2 * integrand_value
        print(f"sigma_in_h: {sigma_in_h}")
        return sigma_in_h     

    def free_way_until_next_interaction(self, mu_T):
        """
        Calculate the distance to the next interaction for the photon.
        
        :param mu_T: Total linear attenuation coefficient at the photon's energy.
        :return: Distance to the next interaction.
        """
        U = random.uniform(0, 1)  # Generate a random number U between 0 and 1
        s = -math.log(U) / mu_T

        return s, U
    
    def move_electron(self, s):
        """
        Move the electron a distance s in its current direction.
        
        :param s: Distance to move the photon.
        :return: A boolean indicating if the photon has left the interaction volume.
        """
        tau = random.uniform(0, s)
        
        position = tuple(np.array(self.position) + tau * np.array(self.direction))

        return position, tau
    
    def angular_deflection(self, s, sigma_el_h):
        def integrand_deflection(sigma_el_h):
            N = self.medium.number_density()
            lambda_el_1 =1/(2 * N * sigma_el_h * self.A_0 * (1 + self.A_0) * (np.log((self.mu_c + self.A_0) / self.A_0) - self.mu_c / (self.mu_c + self.A_0)))
            return lambda_el_1
        chi = np.arccos(np.exp(-s / integrand_deflection(sigma_el_h)))
        phi = 2 * np.pi * random.uniform(0, 1)
        return chi, phi
    
    def calculate_energy_loss(self, s):
        # Constantes y valores del medio
        N = self.medium.number_density()
        Z = self.medium.atomic_number
        I = self.medium.I
        mc2 = 510.999  # Energía de reposo del electrón en keV
        re = 2.81794e-13
        
        # Calculamos los valores de beta y gamma
        gamma = (self.energy + mc2) / mc2
        beta = np.sqrt((gamma**2 - 1) / gamma**2)


        # Corrección de densidad de Fermi 
        delta_f = 0 

        S_s_E = N * Z * 2 * np.pi * re**2 * mc2 / beta**2 * (np.log((self.energy**2/I**2)* (gamma + 1) / 2) + 1 - beta**2 - (2 * gamma - 1)/ gamma**2 * np.log(2) + 1/8 * ((gamma - 1) / gamma)**2 - delta_f)

        # La pérdida de energía w es simplemente S_s(E) multiplicado por la distancia s
        w = S_s_E * s
        
        if w >= self.energy:
            return self.energy
        else:
            return w

            
    def electron_simulation(self):
        continue_simulation = True
        sigma_el_h = self.effective_section_el()
        sigma_in_h = self.effective_section_in()
        mu_t = self.attenuation_coefficient_calculation(sigma_el_h, sigma_in_h)
        print(f"mu_t: {mu_t}")
        s, U = self.free_way_until_next_interaction(mu_t)
        print(f"s: {s}, U: {U}")
        self.position, tau = self.move_electron(s)
        print(f"Position: {self.position}")
        if self.position[2] > 0:
            chi, phi = self.angular_deflection(s, sigma_el_h)
            print(f'Chi:{chi}')
            w = self.calculate_energy_loss(s)
            print(f"Energy lost: {w}")
            self.direction = self.geometry.rotate_vector(self.direction_0, chi, phi)
            print(f"New direction: {self.direction}")
            self.energy = self.energy - w
            if self.energy < self.medium.E_abs_el:
                continue_simulation = False
            else:
                self.positon = self.move_electron(s - tau)
                print(f"New position: {self.position}")
                if self.position[2] < 0:
                    print(f'Electron exited the permitted space at position {self.position}...')
                    continue_simulation = False
                else:
                    
        else:
            print(f'Electron exited the permitted space at position {self.position}...')
            continue_simulation = False
        
