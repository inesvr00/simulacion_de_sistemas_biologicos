import numpy as np


class InteractionVolume:
    def __init__(self, origin=(0, 0, 0), direction=(0, 0, 1)):
        """
        Initialize the InteractionVolume with an origin and a direction.

        Parameters:
        origin (tuple): The starting point of the volume, default is (0, 0, 0).
        direction (tuple): The direction of the volume, default is (0, 0, 1).
        """
        
        self.origin = origin
        self.direction = direction

    def path_length(self, point):
        """
        Calculate the path length of a photon from the origin to the point where
        it exits the volume or is absorbed.

        Parameters:
        point (tuple): The point where the photon exits or is absorbed.

        Returns:
        float: The path length from the origin to the exit/absorption point. 
        Returns 0 if the point is not within the volume.
        """

        if not self.contains(point):
            return 0
        else:
            return point[2] - self.origin[2]  # Photon moves in the z direction
        
    def rotate_vector(self, original_vector, theta, phi):
        """
        Rotate a vector by given angles theta and phi.

        Parameters:
        original_vector (array-like): The original vector to be rotated.
        theta (float): The angle to rotate around the x-axis.
        phi (float): The angle to rotate around the z-axis.

        Returns:
        numpy.ndarray: The new vector after rotation.
        """
        
        original_vector = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0])
        # Extract individual components from the original vector
        u, v, w = original_vector

        # Calculate new vector components after rotation
        u_prime = u * np.cos(theta) + (v * np.cos(phi) - w * np.sin(phi)) * np.sin(theta) / np.sqrt(1 - w**2)
        v_prime = v * np.cos(theta) + (w * np.cos(phi) + u * np.sin(phi)) * np.sin(theta) / np.sqrt(1 - w**2)
        w_prime = w * np.cos(theta) - np.sqrt(1 - w**2) * np.sin(theta) * np.cos(phi)
    
        new_direction = np.array([u_prime, v_prime, w_prime])

        return new_direction
