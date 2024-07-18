import numpy as np
import matplotlib.pyplot as plt
import random

# Defining the theoretical solutions for the parent and daughter nuclei
def theoretical_solution(N_p0, N_h0, lambda_p, lambda_h, t):
    N_p_t = N_p0 * np.exp(-lambda_p * t)
    term1 = N_p0 * (lambda_p / (lambda_h - lambda_p)) * (np.exp(-lambda_p * t) - np.exp(-lambda_h * t))
    term2 = N_h0 * np.exp(-lambda_h * t)
    N_h_t = term1 + term2
    return N_p_t, N_h_t

class RadioactiveDecayChainSimulation:
    def __init__(self, N_p0, N_h0, N_h20, decay_constants, dt=0.1):
        self.N_p0 = N_p0
        self.N_h0 = N_h0
        self.N_h20 = N_h20
        self.lambda_p, self.lambda_h = decay_constants
        self.dt = dt
        self.time = 0
        self.N_p = N_p0
        self.N_h = N_h0
        self.N_h2 = N_h20
        self.N_ht_history = self.N_h + self.N_h2
        self.time_history = [self.time]
        self.N_p_history = [self.N_p]
        self.N_h_history = [self.N_h]
        self.N_h2_history = [self.N_h2]
        self.N_ht_history = [self.N_ht_history]
        self.theoretical_N_p_history = [self.N_p]
        self.theoretical_N_h_history = [self.N_h]

    def decay_step_stochastic(self):
        for _ in range(self.N_p):
            if np.random.rand() < 1 - np.exp(-self.lambda_p * self.dt):
                self.N_p -= 1
                self.N_h += 1
        
        for _ in range(self.N_h):
            if np.random.rand() < 1 - np.exp(-self.lambda_h * self.dt):
                self.N_h -= 1
                self.N_h2 += 1

    def run_simulation(self, total_time):
        steps = int(total_time / self.dt)
        for _ in range(steps):
            self.decay_step_stochastic()
            self.time += self.dt
            self.time_history.append(self.time)
            self.N_p_history.append(self.N_p)
            self.N_h_history.append(self.N_h)
            self.N_h_history.append(self.N_h2)
            # Calculate the theoretical values at this time
            N_p_t, N_h_t = theoretical_solution(self.N_p0, self.N_h0, self.lambda_p, self.lambda_h, self.time)
            self.theoretical_N_p_history.append(N_p_t)
            self.theoretical_N_h_history.append(N_h_t)

    def plot_simulation(self):
        plt.plot(self.time_history, self.N_p_history, label='Simulated N_p (Stochastic)')
        plt.plot(self.time_history, self.N_h_history, label='Simulated N_h (Stochastic)')
        plt.plot(self.time_history, self.N_h2_history, label='Simulated N_h2 (Stochastic)')
        plt.plot(self.time_history, self.theoretical_N_p_history, '--', label='Theoretical N_p')
        plt.plot(self.time_history, self.theoretical_N_h_history, '--', label='Theoretical N_h')
        plt.xlabel('Time')
        plt.ylabel('Number of Nuclei')
        plt.title('Radioactive Decay Chain Stochastic Simulation vs Theory')
        plt.legend()
        plt.grid(True)
        plt.show()

N_p0 = 100000  
N_h0 = 0       
N_h20 = 0
decay_constants = [0.693 / 5, 0.693 / 10]  
total_time = 100  

simulation = RadioactiveDecayChainSimulation(N_p0, N_h0, N_h20, decay_constants)
simulation.run_simulation(total_time)
simulation.plot_simulation()
