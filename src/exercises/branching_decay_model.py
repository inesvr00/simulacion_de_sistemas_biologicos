import numpy as np
import random
import matplotlib.pyplot as plt

def theoretical_solution_branching(N_p, N_a, N_b, lambda_a, lambda_b, dt):
    N_p_t = N_p * (1 -(lambda_a + lambda_b) * dt)
    N_a_t = N_a + lambda_a * N_p * dt
    N_b_t = N_b + lambda_b * N_p * dt
    return N_p_t, N_a_t, N_b_t

class RadioactiveDecayBranchingSimulation:
    def __init__(self, N_p0, decay_constants, dt=0.1):
        self.lambda_a, self.lambda_b = decay_constants
        self.dt = dt
        self.time = 0
        self.time_history = [self.time]
        self.N_p = N_p0
        self.N_a = 0
        self.N_b = 0
        self.N_p_history = [self.N_p]
        self.N_a_history = [self.N_a]
        self.N_b_history = [self.N_b]
        self.theoretical_N_p_history = [self.N_p]
        self.theoretical_N_a_history = [self.N_a]
        self.theoretical_N_b_history = [self.N_b]

    def decay_nuclei_branching(self):
        decayed_to_a = 0
        decayed_to_b = 0
        for _ in range(self.N_p):
            if random.random() < (1 - np.exp(-(self.lambda_a + self.lambda_b) * self.dt)):
                if random.random() < self.lambda_a / (self.lambda_a + self.lambda_b):
                    decayed_to_a += 1
                else:
                    decayed_to_b += 1
        self.N_p -= (decayed_to_a + decayed_to_b)
        self.N_a += decayed_to_a
        self.N_b += decayed_to_b

    def run_simulation_step(self):
        self.decay_nuclei_branching()
        self.time += self.dt
        self.time_history.append(self.time)
        self.N_p_history.append(self.N_p)
        self.N_a_history.append(self.N_a)
        self.N_b_history.append(self.N_b)
        # Calculate the theoretical values at this time
        N_p_t, N_a_t, N_b_t = theoretical_solution_branching(self.N_p, self.N_a, self.N_b, self.lambda_a, self.lambda_b, self.dt)
        self.theoretical_N_p_history.append(N_p_t)
        self.theoretical_N_a_history.append(N_a_t)
        self.theoretical_N_b_history.append(N_b_t)

    def run_simulation(self):
        while np.any(self.N_p > 0):  
            self.run_simulation_step()


    def plot_simulation(self):
        plt.plot(self.time_history, self.N_p_history, label='Núcleos padre')
        plt.plot(self.time_history, self.N_a_history, label='Núcleos hijo A')
        plt.plot(self.time_history, self.N_b_history, label='Núcleos hijo B')
        plt.plot(self.time_history, self.theoretical_N_p_history, '--', label='Núcleos padre teórico')
        plt.plot(self.time_history, self.theoretical_N_a_history, '--', label='Núcleos hijo A teórico')
        plt.plot(self.time_history, self.theoretical_N_b_history, '--', label='Núcleos hijo B teórico')
        plt.xlabel('Tiempo/s')
        plt.ylabel('Número de núcleos')
        plt.title('Variación del estado de los núcleos en el tiempo')
        plt.legend()
        plt.grid(True)
        plt.savefig('/home/ines/simulacion_de_sistemas_biologicos/out/3-ejercicio/radioactive_decay_branching_simulation_distribution.png')
        plt.show()

if __name__ == "__main__":
    N_p0 = 10000
    decay_constants = [0.693 / 5, 0.693 / 10]
    simulation = RadioactiveDecayBranchingSimulation(N_p0, decay_constants)
    simulation.run_simulation()
    simulation.plot_simulation()