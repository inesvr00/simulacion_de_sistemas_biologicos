import numpy as np


class InteractionVolume:
    def __init__(self, origin=(0, 0, 0), direction=(0, 0, 1)):
        self.origin = origin
        self.direction = direction

    def contains(self, point):
        """
        Check if the given point is within the interaction volume.
        The volume is defined as all points with z > 0.
        """
        x, y, z = point
        return z > 0

    def path_length(self, point):
        """
        Simulate the path of a photon until it is absorbed or exits the volume.
        This function returns the path length from the origin to the point where
        the photon exits the volume or is absorbed, assuming it moves in a straight
        line in the z-direction.
        """
        if not self.contains(point):
            return 0
        else:
            return point[2] - self.origin[2]  # Assuming photon moves in the z direction
        
    def rotate_vector(self, original_vector, theta, phi):
        original_vector = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0])
        # Extraer las componentes individuales del vector original.
        u, v, w = original_vector

        # Cálculo de los nuevos componentes del vector después de la rotación.
        u_prime = u * np.cos(theta) + (v * np.cos(phi) - w * np.sin(phi)) * np.sin(theta) / np.sqrt(1 - w**2)
        v_prime = v * np.cos(theta) + (w * np.cos(phi) + u * np.sin(phi)) * np.sin(theta) / np.sqrt(1 - w**2)
        w_prime = w * np.cos(theta) - np.sqrt(1 - w**2) * np.sin(theta) * np.cos(phi)
    
        # La nueva dirección del vector como array de numpy.
        new_direction = np.array([u_prime, v_prime, w_prime])

        return new_direction
