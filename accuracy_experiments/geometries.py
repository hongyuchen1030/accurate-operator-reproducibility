import numpy as np

# Geometric functions
class angle:
    # Defines an angle class to allow easy switching between degrees and radians
    def __init__(self, theta: float, radians=True):
        if radians:
            self.rad = theta
            self.deg = theta * 180 / np.pi
        else:
            self.deg = theta
            self.rad = theta * np.pi / 180

    def __add__(self, other):
        return angle(self.rad + other.rad)

    def __sub__(self, other):
        return angle(self.rad - other.rad)


class coordinates:
    # The coordinates class allows for conversion between lat-lon and xyz (used for calculations)
    def __init__(self, x, R=1, degrees=True, lonlatorder=True):
        if len(x) == 2:  # lat-lon inputs
            if lonlatorder:
                self.lon = angle(x[0], radians=not degrees)
                self.lat = angle(x[1], radians=not degrees)
            else:
                self.lon = angle(x[1], radians=not degrees)
                self.lat = angle(x[0], radians=not degrees)
            self.R = R
            self.latlon2xyz()
        elif len(x) == 3:  # xyz inputs
            self.x = x
            self.R = np.linalg.norm(x)
            self.xyz2latlon()
        else:
            raise ValueError(f"ERROR: Cannot evaluate an x coordinate of length {len(x)} in coordinate class")

    def latlon2xyz(self):
        # Converts lat-lon to xyz coordinates on a sphere of radius R
        self.x = self.R * np.array([
            np.cos(self.lat.rad) * np.cos(self.lon.rad),
            np.cos(self.lat.rad) * np.sin(self.lon.rad),
            np.sin(self.lat.rad)
        ])

    def xyz2latlon(self):
        # Converts xyz coordinates to latitude-longitude
        self.lat = angle(np.arctan2(self.x[2], np.sqrt(self.x[0] ** 2 + self.x[1] ** 2)))

        self.lon = angle(np.arctan2(self.x[1], self.x[0]))

    def latlon2xy(self, origin: "coordinates") -> np.ndarray:
        # Returns the location of the point on planar projection relative to the origin
        dlat = self.lat - origin.lat
        dlon = self.lon - origin.lon
        dx = dlon.rad * self.R * np.cos(origin.lat.rad)
        dy = dlat.rad * self.R
        return np.array([dx, dy])


def rot_mat(theta: float, ax: int) -> np.ndarray:
    # Returns the rotation matrix when rotating by an angle theta around axis ax (0, 1, 2 for x, y, z)
    c = np.cos(theta)
    s = np.sin(theta)
    if ax == 0:  # Rotation around x-axis
        M = np.array([[1, 0,  0],
                      [0, c, -s],
                      [0, s,  c]])
    elif ax == 1:  # Rotation around y-axis
        M = np.array([[ c, 0, s],
                      [ 0, 1, 0],
                      [-s, 0, c]])
    elif ax == 2:  # Rotation around z-axis
        M = np.array([[c, -s, 0],
                      [s,  c, 0],
                      [0,  0, 1]])
    else:
        raise ValueError("Axis must be 0, 1, or 2 for x, y, or z rotation")
    return M


# Example usage or setup of psi for further computations
psi = np.linspace(0, 2 * np.pi, 101)
