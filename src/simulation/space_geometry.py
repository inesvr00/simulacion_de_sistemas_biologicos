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
        """
        Rotate the original vector by the given angles.
        
        :param original_vector: The original vector to rotate.
        :param theta: The polar angle.
        :param phi: The azimuthal angle.
        :return: The rotated vector.
        """
        # Convert from spherical to Cartesian coordinates
        r = np.linalg.norm(original_vector)
        stheta = r * np.cos(theta)
        r_xy = r * np.sin(theta)
        x = r_xy * np.cos(phi)
        y = r_xy * np.sin(phi)
        z = stheta

        # The new Cartesian coordinates are the rotated vector
        rotated_vector = np.array([x, y, z])
        return rotated_vector
