import numpy as np
import random
import matplotlib.pyplot as plt

def theoretical_solution(N_p, N_h1, lambda_p, lambda_h, dt):
    N_p_t = N_p * (1 - lambda_p * dt)
    N_h_t = N_h1 + (lambda_p * N_p * dt) - (lambda_h * N_h1 * dt)
    return N_p_t, N_h_t

class RadioactiveDecayChainSimulation:
    def __init__(self, num_nuclei, decay_constants):
        self.num_nuclei = num_nuclei
        self.N_p0 = num_nuclei
        self.N_p = self.N_p0
        self.N_h1 = 0
        self.N_h2 = 0
        self.lambda_p, self.lambda_h = decay_constants
        
        # Initial state of nuclei: 2 for parent, 1 for first daughter, 0 for second daughter
        self.state = np.ones(num_nuclei, dtype=int) * 2
        self.time = 0
        
        # Histories for plotting and analysis
        self.time_history = [self.time]
        self.remaining_nuclei_history = [self.num_nuclei]
        self.state_distribution_history = [[np.sum(self.state == i) for i in range(0, 3)]]
        self.theoretical_N_p_history = [self.N_p]
        self.theoretical_N_h_history = [self.N_h1]


    def decay_nuclei_chain(self, dt):
        for i in range(self.num_nuclei):
            current_state = self.state[i]
            if current_state == 2:
                if random.random() < (1 - np.exp(-self.lambda_p * dt)):
                    self.state[i] = 1
                    self.N_p -= 1
                    self.N_h1 += 1
            elif current_state == 1:
                if random.random() < (1 - np.exp(-self.lambda_h * dt)):
                    self.state[i] = 0
                    self.N_h2 += 1
                    self.N_h1 -= 1
 
    def run_simulation_step(self):
        dt = 0.1
        self.decay_nuclei_chain(dt)
        self.time += dt  
        self.time_history.append(self.time)

        # Update histories
        self.remaining_nuclei_history.append(np.count_nonzero(self.state > 0))
        self.state_distribution_history.append([np.sum(self.state == i) for i in range(0, 3)])
        
        # Update theoretical histories
        N_p_t, N_h_t = theoretical_solution(self.N_p, self.N_h1, self.lambda_p, self.lambda_h, dt)
        self.theoretical_N_p_history.append(N_p_t)
        self.theoretical_N_h_history.append(N_h_t)

    def run_simulation(self):
        while np.any(self.state >= 1):  
            self.run_simulation_step()

    def plot_simulation(self):
        plt.plot(self.time_history, self.remaining_nuclei_history)
        plt.xlabel('Time')
        plt.ylabel('Remaining Initial Radioactive Nuclei')
        plt.title('Radioactive Decay Chain Simulation')
        plt.grid(True)
        plt.savefig('/home/ines/simulacion_de_sistemas_biologicos/out/3-ejercicio/radioactive_decay_chain_simulation.png')
        plt.show()
        
    def plot_state_distributions(self):
        plt.figure(figsize=(10, 6))
        for state_index, state_values in enumerate(zip(*self.state_distribution_history)):
            plt.plot(self.time_history, state_values, label=f'Estado {state_index}')
        plt.plot(self.time_history, self.theoretical_N_p_history, '--', label='Teórica estado 2')
        plt.plot(self.time_history, self.theoretical_N_h_history, '--', label='Teórica estado 1')
        plt.xlabel('Tiempo/s')
        plt.ylabel('Número de núcleos')
        plt.title('Distribución de estado de los núcleos en el tiempo')
        plt.legend()
        plt.grid(True)
        plt.savefig('/home/ines/simulacion_de_sistemas_biologicos/out/3-ejercicio/radioactive_decay_chain_simulation_distribution.png')
        plt.show()

if __name__ == "__main__":
    num_nuclei = 10000
    decay_constants = [0.693 / 5, 0.693 / 10]
    simulation = RadioactiveDecayChainSimulation(num_nuclei, decay_constants)
    simulation.run_simulation()
    simulation.plot_simulation()
    simulation.plot_state_distributions()
