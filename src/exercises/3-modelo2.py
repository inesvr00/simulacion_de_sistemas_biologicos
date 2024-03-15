import numpy as np
import random
import matplotlib.pyplot as plt

class RadioactiveDecaySimulation:
    def __init__(self, num_nuclei, decay_constant):
        self.num_nuclei = num_nuclei
        self.decay_constant = decay_constant
        self.state = np.ones(num_nuclei, dtype=int)
        self.time = 0
        self.time_history = []
        self.num_nuclei_history = []
        
    def time_step_determination(self):
        rate = self.decay_constant * np.sum(self.state)
        if rate > 0:
            random_number = np.random.uniform(low=1e-10, high=1.0)
            time_step = -np.log(random_number) / rate
            return time_step
        else:
            return None  # Or a suitable alternative to indicate the simulation should end
        
    def calculate_dt(self):
        return -np.log(1 - random.random) / self.decay_constant[0]
    
    def decay_nucleus(self):
        return 0 
    
    def run_simulation_step(self):
        N_i = np.sum(self.state)
        if N_i == 0:
            return 0 
        random_index = np.random.choice(np.where(self.state == 1)[0])
        new_state = np.copy(self.state)
        new_state[random_index] = self.decay_nucleus()
        self.state = new_state
        self.time_step = self.time_step_determination()
        if self.time_step is not None:  # This check is now somewhat redundant since we're returning 0 instead of None, but it's good practice
            self.time += self.time_step
        self.time_history.append(self.time)
        self.num_nuclei_history.append(N_i)
        return N_i

    def run_simulation(self):
        total_steps = 0
        while True:
            N_i = self.run_simulation_step()
            if N_i == 0:
                break
            total_steps += 1
            print(f"Step {total_steps}: Remaining nuclei: {N_i}, Time: {self.time:.2f} s")
            
    def plot_simulation(self):
        theoretical_time = np.linspace(0, self.time, len(self.num_nuclei_history))
        theoretical_remaining_nuclei = self.num_nuclei * np.exp(-self.decay_constant * theoretical_time)
        plt.plot(self.time_history, self.num_nuclei_history, label='Simulation')  # Corrected this line
        plt.plot(theoretical_time, theoretical_remaining_nuclei, label='Theoretical', linestyle='--')
        plt.xlabel('Time')
        plt.ylabel('Remaining Nuclei')
        plt.title('Radioactive Decay Simulation vs Theoretical Curve')
        plt.legend()
        plt.grid(True)
        plt.savefig('radioactive_decay_simulation.png')
        plt.show()

        
if __name__ == "__main__":
    num_nuclei = 10000
    decay_constant = 0.1315
    simulation = RadioactiveDecaySimulation(num_nuclei, decay_constant)
    simulation.run_simulation()
    simulation.plot_simulation()                
        