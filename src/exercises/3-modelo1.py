import numpy as np
import random
from loguru import logger
import matplotlib.pyplot as plt

class RadioactiveDecaySimulation:
    def __init__(self, num_nuclei, decay_constant):
        self.num_nuclei = num_nuclei
        self.decay_constant = decay_constant
        self.state = np.ones(num_nuclei, dtype=int)
        self.time = 0
        self.remaining_nuclei_history = []

    def decay_probability(self, dt):
        return 1 - np.exp(-self.decay_constant * dt)

    def calculate_dt(self, p):
        return -np.log(1 - p) / self.decay_constant

    def decay_nucleus(self, nucleus):
        return 0 if random.random() <= nucleus else 1

    def run_simulation_step(self):
        N_i = np.sum(self.state)
        p = self.decay_probability(self.calculate_dt(self.decay_constant))
        new_state = np.copy(self.state)
        for i in range(self.num_nuclei):
            if self.state[i] == 1:
                new_state[i] = self.decay_nucleus(p)
        self.state = new_state
        self.time += self.calculate_dt(p)
        self.remaining_nuclei_history.append(N_i)
        return N_i

    def run_simulation(self):
        total_steps = 0
        while True:
            N_i = self.run_simulation_step()
            if N_i == 0:
                break
            total_steps += 1
            logger.info(f"Step {total_steps}: Remaining nuclei: {N_i}, Time: {self.time:.2f} s")
            
    def plot_simulation(self):
        theoretical_time = np.linspace(0, self.time, len(self.remaining_nuclei_history))
        theoretical_remaining_nuclei = self.num_nuclei * np.exp(-self.decay_constant * theoretical_time)
        plt.plot(np.arange(len(self.remaining_nuclei_history)), self.remaining_nuclei_history, label='Simulation')
        plt.plot(theoretical_time, theoretical_remaining_nuclei, label='Theoretical', linestyle='--')
        plt.xlabel('Time Step')
        plt.ylabel('Remaining Nuclei')
        plt.title('Radioactive Decay Simulation vs Theoretical Curve')
        plt.legend()
        plt.grid(True)
        plt.show()

if __name__ == "__main__":
    num_nuclei = 100000
    decay_constant = 0.1315
    simulation = RadioactiveDecaySimulation(num_nuclei, decay_constant)
    simulation.run_simulation()
    simulation.plot_simulation()
