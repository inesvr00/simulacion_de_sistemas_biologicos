import numpy as np
import random
import matplotlib.pyplot as plt

class RadioactiveDecayChainSimulation:
    def __init__(self, num_nuclei, decay_constants, mode):
        self.num_nuclei = num_nuclei
        self.decay_constants = decay_constants
        self.state = np.ones(num_nuclei, dtype=int) * 2
        self.mode = mode
        self.time = 0
        self.time_history = [self.time]
        self.time_history = [self.time]
        self.remaining_nuclei_history = [self.num_nuclei]
        self.state_distribution_history = [[np.sum(self.state == i) for i in range(1, 5)]]

    def decay_nuclei_chain(self, dt):
        for i in range(self.num_nuclei):
            current_state = self.state[i]
            if current_state == 2:
                if random.random() < (1 - np.exp(-self.decay_constants[0] * dt)):
                    self.state[i] = 1
            elif current_state == 1:
                if random.random() < (1 - np.exp(-self.decay_constants[1] * dt)):
                    self.state[i] = 0

            
    def decay_nuclei_branching(self, dt):
        for i in range(self.num_nuclei):
            current_state = self.state[i]
            if current_state == 2:
                # State 4 can decay to state 3 or to state 1
                if random.random() < (1 - np.exp(-self.decay_constants[0] * dt)):
                    self.state[i] = 1
                    self.state[i] = 0
    
    def run_simulation_step(self):
        dt = 0.1
        if self.mode == 1:
            self.decay_nuclei_chain(dt)
        elif self.mode == 2:
            self.decay_nuclei_branching(dt)
        self.time += dt  
        self.time_history.append(self.time)
        self.remaining_nuclei_history.append(np.count_nonzero(self.state > 0))
        self.state_distribution_history.append([np.sum(self.state == i) for i in range(0, 3)])

    def run_simulation(self):
        total_steps = 0
        while np.any(self.state >= 1):  
            self.run_simulation_step()
            total_steps += 1

    def plot_simulation(self):
        plt.plot(self.time_history, self.remaining_nuclei_history)
        plt.xlabel('Time')
        plt.ylabel('Remaining Initial Radioactive Nuclei')
        plt.title('Radioactive Decay Chain Simulation')
        plt.grid(True)
        plt.savefig('radioactive_decay_chain_simulation.png')
        plt.show()
        
    def plot_state_distributions(self):
        plt.figure(figsize=(10, 6))
        for state_index, state_values in enumerate(zip(*self.state_distribution_history)):
            plt.plot(self.time_history, state_values, label=f'State {state_index}')
        plt.xlabel('Time')
        plt.ylabel('Number of Nuclei')
        plt.title('Nuclei State Distribution Over Time')
        plt.legend()
        plt.grid(True)
        plt.savefig('radioactive_decay_chain_simulation_distribution.png')
        plt.show()

if __name__ == "__main__":
    num_nuclei = 100000
    # decay_constants = [0.693 / 20.8, 0.693 / 10]
    decay_constants = [0.693 / 5, 0.693 / 10]
    mode = int(input("Choose a mode: \n 1. Chain Decay \n 2. Branching decay \n" ))
    simulation = RadioactiveDecayChainSimulation(num_nuclei, decay_constants, mode)
    simulation.run_simulation()
    simulation.plot_simulation()
    simulation.plot_state_distributions()
