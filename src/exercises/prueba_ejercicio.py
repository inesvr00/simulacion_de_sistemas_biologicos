import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

class RadioactiveDecay:
    def __init__(self, N_p0, N_h0, N_a0, N_b0, lambda_p, lambda_a, t_range):
        self.N_p0 = N_p0
        self.N_h0 = N_h0
        self.N_a0 = N_a0
        self.N_b0 = N_b0
        self.lambda_p = lambda_p
        self.lambda_a = lambda_a
        self.t_range = t_range

    def decay_simple(self, y, t, lambda_p):
        N_p, N_h = y
        dNp_dt = -lambda_p * N_p
        dNh_dt = lambda_p * N_p
        return [dNp_dt, dNh_dt]

    def decay_branching(self, y, t, lambda_p, lambda_a):
        N_p, N_a, N_b = y
        dNp_dt = -(lambda_p + lambda_a) * N_p
        dNa_dt = lambda_p * N_p
        dNb_dt = lambda_a * N_p
        return [dNp_dt, dNa_dt, dNb_dt]

    def solve_simple_decay(self):
        solution = odeint(self.decay_simple, [self.N_p0, self.N_h0], self.t_range, args=(self.lambda_p,))
        return solution

    def solve_branching_decay(self):
        solution = odeint(self.decay_branching, [self.N_p0, self.N_a0, self.N_b0], self.t_range, args=(self.lambda_p, self.lambda_a))
        return solution

    def plot_results(self, solution_simple, solution_branching, t_range):
        # Plotting the results
        plt.figure(figsize=(12, 5))

        # Desintegración simple
        plt.subplot(1, 2, 1)
        plt.plot(t_range, solution_simple[:, 0], label='Núcleo Padre (N_p)')
        plt.plot(t_range, solution_simple[:, 1], label='Núcleo Hijo (N_h)')
        plt.title('Desintegración Simple')
        plt.xlabel('Tiempo')
        plt.ylabel('Número de Núcleos')
        plt.legend()

        # Desintegración ramificada
        plt.subplot(1, 2, 2)
        plt.plot(t_range, solution_branching[:, 0], label='Núcleo Padre (N_p)')
        plt.plot(t_range, solution_branching[:, 1], label='Núcleo Hijo A (N_a)')
        plt.plot(t_range, solution_branching[:, 2], label='Núcleo Hijo B (N_b)')
        plt.title('Desintegración Ramificada')
        plt.xlabel('Tiempo')
        plt.ylabel('Número de Núcleos')
        plt.legend()

        plt.tight_layout()
        plt.savefig('decay_plots.png')
        plt.close()

# Values and constants
N_p0 = 2000.0
N_h0 = 0.0
N_a0 = 0.0
N_b0 = 0.0
lambda_p = 1
lambda_a = 0.5
t = np.linspace(0, 10, 100)

# Create instance of the class
decay_simulation = RadioactiveDecay(N_p0, N_h0, N_a0, N_b0, lambda_p, lambda_a, t)

# Solve the equations
sol_simple = decay_simulation.solve_simple_decay()
sol_branching = decay_simulation.solve_branching_decay()

# Plot and save the results
decay_simulation.plot_results(sol_simple, sol_branching, t)


