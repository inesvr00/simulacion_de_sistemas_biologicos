from matplotlib import pyplot as plt
import numpy as np

class ThetaDistribution:
    def __init__(self):
        pass
    
    def pdf(self, theta):
        return np.sin(theta) / 2

class PhiDistribution:
    def __init__(self):
        pass
    
    def pdf(self, phi):
        return np.ones_like(phi) / (2 * np.pi)

class DistributionPlotter:
    def __init__(self, num_samples):
        self.num_samples = num_samples
        
    def plot_theta_distribution(self):
        random_sample_u = np.random.uniform(0, 1, self.num_samples)
        theta_samples = np.arccos(1 - 2 * random_sample_u)

        plt.hist(theta_samples, bins=50, density=True, alpha=0.7, color='blue')
        theta_range = np.linspace(0, np.pi, 100)
        plt.plot(theta_range, ThetaDistribution().pdf(theta_range), 'r--', linewidth=2)
        plt.title('Distribución de $\\Theta$')
        plt.xlabel('$\\Theta$ (radianes)')
        plt.ylabel('Densidad de probabilidad')
        plt.savefig('histogram_path_distribution_theta.png')
        plt.show()
    
    def plot_phi_distribution(self):
        phi_samples = np.random.uniform(0, 2 * np.pi, self.num_samples)

        plt.hist(phi_samples, bins=50, density=True, alpha=0.7, color='green')
        phi_range = np.linspace(0, 2 * np.pi, 100)
        plt.plot(phi_range, PhiDistribution().pdf(phi_range), 'r--', linewidth=2)
        plt.title('Distribución de $\\Phi$')
        plt.xlabel('$\\Phi$ (radianes)')
        plt.ylabel('Densidad de probabilidad')
        plt.savefig('histogram_path_distribution_phi.png')
        plt.show()
        
    def plot_random_emission_direction(self):
        # Muestreo de un único par de ángulos
        random_sample_u = np.random.uniform(0, 1)
        theta = np.arccos(1 - 2 * random_sample_u)
        phi = np.random.uniform(0, 2 * np.pi)
        
        # Conversión a coordenadas cartesianas
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        
        # Crear figura y eje 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Dibujar el vector de emisión
        ax.quiver(0, 0, 0, x, y, z, length=1.0, color='b', arrow_length_ratio=0.1)
        
        # Establecer límites y etiquetas
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Título y visualización
        ax.set_title('Dirección aleatoria de emisión en 3D')
        plt.savefig('random_emission_direction.png')
        plt.show()

def main():
    num_samples = 10000
    distribution_plotter = DistributionPlotter(num_samples)
    distribution_plotter.plot_theta_distribution()
    plt.close()
    distribution_plotter.plot_phi_distribution()
    plt.close()
    distribution_plotter.plot_random_emission_direction()

if __name__ == "__main__":
    main()
