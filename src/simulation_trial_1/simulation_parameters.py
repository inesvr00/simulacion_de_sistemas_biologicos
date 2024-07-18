class SimulationParameters:
    def __init__(self, theta_c=20, W_ci=1):
        """
        Initialize the SimulationParameters with the properties of the simulation.

        :param theta_c: Cut angle for elastic collisions.
        :param W_ci: Cut energy for inelastic collisions.
        """
        self.theta_c = theta_c
        self.W_ci = W_ci