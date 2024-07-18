import math
import random
import numpy as np
import pandas as pd
import scipy
from src.simulation_trial_1.comprobacion_perdida_energia_paso_suave import añadir_perdida_energia_paso_suave
from src.simulation_trial_1.kerma_dosis import add_dosis
from src.simulation_trial_1.medium import AbsorbentMedium
from src.simulation_trial_1.space_geometry import InteractionVolume



class ElectronProperties:
    def __init__(self, E0, position_0, direction_0, theta_c=20, W_ci=1, medium = AbsorbentMedium(), geometry = InteractionVolume()):
        """
        Initialize the ElectronProperties with the properties of an electron.

        :param E0: Initial energy of the electron.
        :param position_0: Initial position of the electron.
        :param direction_0: Initial direction of the electron.
        :param theta_c: Cut angle for elastic collisions, default is 20 degrees.
        :param W_ci: Cut energy for inelastic collisions, default is 1.
        :param medium: The medium in which the electron is traveling, default is AbsorbentMedium.
        :param geometry: The geometry of the interaction volume, default is InteractionVolume.
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
        self.N =  self.medium.number_density()
        self.mc2 = 510.999
        
    def __str__(self):
        return f"Electron(Energy={self.energy:.2f} keV, Position={self.position}, Direction={self.direction})"

    def attenuation_coefficient_calculation(self, sigma_el_h, sigma_in_h ):
        """
        Calculate the total attenuation coefficient.

        :param sigma_el_h: Effective cross-section for elastic scattering.
        :param sigma_in_h: Effective cross-section for inelastic scattering.
        :return: Total attenuation coefficient.
        """
        
        mu_T = self.N * (sigma_el_h + sigma_in_h)
        return mu_T
        
    def effective_section_el(self):
        """
        Calculate the effective cross-section for elastic scattering.

        :return: Effective cross-section for elastic scattering.
        """

        sigma_el = (1.290e-14) / (1000 * self.energy)**0.730
        sigma_el_h = sigma_el * self.A_0 * (1 - self.mu_c) / (self.mu_c + self.A_0)
        # print(f"sigma_el_h: {sigma_el_h}")
        return sigma_el_h
    
    def effective_section_in(self):
        """
        Calculate the effective cross-section for inelastic scattering.

        :return: Effective cross-section for inelastic scattering.
        """

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
        # print(f"sigma_in_h: {sigma_in_h}")
        return sigma_in_h     

    def free_way_until_next_interaction(self, mu_T):
        """
        Calculate the distance to the next interaction for the photon.
        
        :param mu_T: Total linear attenuation coefficient at the photon's energy.
        :return: Distance to the next interaction.
        """
        U = random.uniform(0, 1)
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
        """
        Calculate the angular deflection of the electron after a collision.

        :param s: Distance traveled by the electron.
        :param sigma_el_h: Effective cross-section for elastic scattering.
        :return: The deflection angles chi and phi.
        """

        def integrand_deflection(sigma_el_h):
            lambda_el_1 =1/(2 * self.N * sigma_el_h * self.A_0 * (1 + self.A_0) * (np.log((self.mu_c + self.A_0) / self.A_0) - self.mu_c / (self.mu_c + self.A_0)))
            return lambda_el_1
        chi = np.arccos(np.exp(-s / integrand_deflection(sigma_el_h)))
        phi = 2 * np.pi * random.uniform(0, 1)
        return chi, phi
    
    def calculate_energy_loss(self, s):
        """
        Calculate the energy loss of the electron after traveling a distance s.

        :param s: Distance traveled by the electron.
        :return: The energy lost by the electron.
        """

        Z = self.medium.atomic_number
        I = self.medium.I
        mc2 = 510.999
        re = 2.81794e-13
        
        gamma = (self.energy + mc2) / mc2
        beta = np.sqrt((gamma**2 - 1) / gamma**2)

        # Corrección de densidad de Fermi 
        delta_f = 0 

        S_s_E = self.N * Z * 2 * np.pi * re**2 * mc2 / beta**2 * (np.log((self.energy**2/I**2)* (gamma + 1) / 2) + 1 - beta**2 - (2 * gamma - 1)/ gamma**2 * np.log(2) + 1/8 * ((gamma - 1) / gamma)**2 - delta_f)

        w = S_s_E * s
        
        if w >= self.energy:
            return self.energy
        else:
            return w
        
    def calculate_p_el(self, sigma_el_h, mu_t):
        """
        Calculate the probability of an elastic collision.

        :param sigma_el_h: Effective cross-section for elastic scattering.
        :param mu_t: Total attenuation coefficient.
        :return: Probability of an elastic collision.
        """

        return self.N * sigma_el_h / mu_t
    
    def calculate_p_in(self, sigma_in_h, mu_t):
        """
        Calculate the probability of an inelastic collision.

        :param sigma_in_h: Effective cross-section for inelastic scattering.
        :param mu_t: Total attenuation coefficient.
        :return: Probability of an inelastic collision.
        """

        return self.N * sigma_in_h / mu_t
    
    def calculate_theta(self, U):
        """
        Calculate the scattering angle theta for an elastic collision.

        :param U: Random variable.
        :return: Scattering angle theta.
        """

        return self.A_0 * U / (1 + self.A_0 - U)

    def calculate_phi_k(self, k_c, k, a):
        """
        Calculate the function phi_k used in inelastic collision calculations.

        :param k_c: Cut energy ratio.
        :param k: Energy ratio.
        :param a: Parameter a.
        :return: Value of phi_k.
        """

        phi_k = (k**-2 + 5 * a) * np.heaviside(k - k_c, 1) * np.heaviside(1/2 - k, 1)
        return phi_k
    
    def calculate_f_in_k(self, k_c, k, a):
        """
        Calculate the function f_in_k used in inelastic collision calculations.

        :param k_c: Cut energy ratio.
        :param k: Energy ratio.
        :param a: Parameter a.
        :return: Value of f_in_k.
        """
        
        f_in_k =(1 / k**2 + 1 / (1 - k)**2 - 1 / (k * (1-k)) + a * (1 + 1 / (k * (1 - k)))) * np.heaviside(k - k_c, 1) * np.heaviside(1/2 - k, 1)
        return f_in_k
    
    def calculate_inellastic_energy_lost(self):
        """
        Calculate the energy lost during an inelastic collision.

        :return: Energy lost and scattering angle theta.
        """

        k_c = self.W_ci / self.energy
        a = (self.energy / (self.energy + self.mc2))**2
        continue_function = True
        p_1 = k_c / (1 - 2 * k_c) * k_c**-2
        p_2 = 2 / (1 - k_c)
        
        while continue_function:
            U_1 = random.uniform(0, 1)
            U_2 = random.uniform(0, 1)

            if U_1 < p_1:
                k = k_c / (1 - U_2 * (1 - 2 * k_c))
            else:
                k = k_c + (U_2 * (1 - 2 * k_c) / 2)

            U_3 = random.uniform(0, 1)

            if U_3 * self.calculate_phi_k(k_c, k, a) > self.calculate_f_in_k(k_c, k, a):
                continue_function = False

        W = k * self.energy
        theta = np.arccos(np.sqrt((self.energy - W) / self.energy * (self.energy + 2 * self.mc2) / (self.energy - W + 2 * self.mc2)))

        return W, theta



    def electron_simulation(self):
        continue_simulation = True
        df_energia_suave = pd.DataFrame(columns=['r', 'energy'])
        df_dosis =  pd.DataFrame(columns=['energy', 'z'])
        n_elastic = 0
        n_inelastic = 0
        while continue_simulation:
            # print(f"Energy electron: {self.energy}")
            sigma_el_h = self.effective_section_el()
            sigma_in_h = self.effective_section_in()
            mu_t = self.attenuation_coefficient_calculation(sigma_el_h, sigma_in_h)
            # print(f"mu_t: {mu_t}")
            s, U = self.free_way_until_next_interaction(mu_t)
            self.position, tau = self.move_electron(s)
            # print(f"Position: {self.position}")
            # if self.position[2] > 0: Descomentar para ejecutar la Simulación 1
            if self.position[2] > 0 and self.position[2] < 0.5: # Descomentar para ejecutar la Simulación 2
                chi, phi = self.angular_deflection(s, sigma_el_h)
                # print(f'Chi:{chi}')
                w = self.calculate_energy_loss(s)
                
                df_energia_suave = añadir_perdida_energia_paso_suave(df_energia_suave, w, self.energy)
                # print(f"Energy lost: {w}")
                self.direction = self.geometry.rotate_vector(self.direction_0, chi, phi)
                df_dosis = add_dosis(df_dosis, w, self.position[2])
                # print(f"New direction: {self.direction}")
                self.energy = self.energy - w
                if self.energy < self.medium.E_abs_el:
                    df_dosis = add_dosis(df_dosis, self.energy, self.position[2])
                    continue_simulation = False
                else:
                    self.positon = self.move_electron(s - tau)
                    # print(f"New position: {self.position}")
                    # if self.position[2] < 0: # Descomentar para ejecutar la Simulación 1
                    if self.position[2] > 0 and self.position[2] < 0.5: # Descomentar para ejecutar la Simulación 2
                        # print(f'Electron exited the permitted space at position {self.position}...')
                        self.energy = 0
                        continue_simulation = False
                    else:
                        p_el = self.calculate_p_el(sigma_in_h, mu_t)
                        p_in = self.calculate_p_in(sigma_in_h, mu_t)
                        if U <= p_el: #Elástico
                            n_elastic += 1
                            theta = self.calculate_theta(U)
                            phi = 2 * np.pi * U
                            self.direction = self.geometry.rotate_vector(self.direction, theta, phi)
                        else: #Inelástica
                            n_inelastic += 1
                            W, theta = self.calculate_inellastic_energy_lost()
                            self.direction = self.geometry.rotate_vector(self.direction, theta, phi)
                            if self.energy > W:
                                df_dosis = add_dosis(df_dosis, W, self.position[2])
                                self.energy = self.energy - W
                            else:
                                self.energy = 0
                                continue_simulation = False
                        if self.energy <= self.medium.E_abs_el:
                            df_dosis = add_dosis(df_dosis, self.energy, self.position[2])
                            continue_simulation = False
            else:
                # print(f'Electron exited the permitted space at position {self.position}...')
                continue_simulation = False
                self.energy = 0
        return df_energia_suave, df_dosis, n_elastic, n_inelastic
        
