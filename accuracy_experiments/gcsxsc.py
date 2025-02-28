import numpy as np
import warnings
from geometries import angle, coordinates
import os
from utils import DecomposedFloat, format_latitude, read_great_circle_arcs_from_csv_exponent, write_great_circle_arcs_to_csv_exponent, format_offset
import csv
class gc_segment:
    # Defines a great circle segment from p_gc_1 to p_gc_2.
    def __init__(self, p_gc_1, p_gc_2, long_segment=False, R_tol=10**-6):
        self.p_1 = p_gc_1
        self.p_2 = p_gc_2
        # Check that the norms of the key vectors all lie on the sphere
        self.R = np.linalg.norm(self.p_1)
        if np.abs(self.R - np.linalg.norm(self.p_2)) / self.R > R_tol:
            warnings.warn("Difference between norms of p_1 and p_2 above tolerance in gc_segment.")
        n_gc = np.cross(self.p_1, self.p_2)
        self.n_hat = n_gc / np.linalg.norm(n_gc)
        if long_segment:
            # Reverse direction of the normal vector for opposite direction
            self.n_hat = -self.n_hat
        # Normalize for ongoing calculations
        self.p_1 = self.p_1 / self.R
        self.p_2 = self.p_2 / self.R
        self.p_1_perp = np.cross(self.n_hat, self.p_1)  # Arc from p_1 to p_2
        # self.beta_2 = self.find_beta(self.p_2)

    def find_beta(self, p2):
        # Finds the angle along the great circle segment from p_1 to p2.
        d1 = np.dot(p2, self.p_1)
        d2 = np.dot(p2, self.p_1_perp)
        if d2 > 0:
            return angle(np.arccos(d1))
        else:
            return angle(2 * np.pi - np.arccos(d1))


class sc:
    def __init__(self, z_0):
        self.n_hat = np.array([0.0, 0.0, 1.0])
        self.r = self.alpha = np.arccos(z_0)
        self.p = z_0


class intersection:
    def __init__(self, gc_segment, sc):
        self.gc_segment = gc_segment
        self.sc = sc

    def N_intersection(self):
        if not hasattr(self, 'N'):
            psi = np.abs(np.pi / 2 - np.arccos(np.dot(self.gc_segment.n_hat, self.sc.n_hat)))
            if psi < self.sc.r:
                self.N = 2
            elif psi == self.sc.r:
                self.N = 1
            else:
                self.N = 0
        return self.N

    def calculate_intersection_line(self):
        # Finds the line of intersection between a great circle and small circle on a sphere.
        nG = self.gc_segment.n_hat
        nS = self.sc.n_hat
        M = np.array([[2,0,0,nG[0],nS[0]],
                        [0,2,0,nG[1],nS[1]],
                        [0,0,2,nG[2],nS[2]],
                        [nG[0],nG[1],nG[2],0,0],
                        [nS[0],nS[1],nS[2],0,0]
        ])
        b = np.array([0, 0, 0, 0, self.sc.p]).transpose()
        self.p_int = np.linalg.lstsq(M, b, rcond=None)[0][:3]
        s = np.cross(self.gc_segment.n_hat, self.sc.n_hat)
        self.s_int_hat = s / np.linalg.norm(s)

    def gc_sc_intersect(self):
        # Calculates points of intersection between the line l = p + s * t and the great circle.
        if not hasattr(self, 'p_int'):
            self.calculate_intersection_line()
        norm_p = np.linalg.norm(self.p_int)
        if norm_p > 1:
            roots = []
        elif norm_p == 1:
            roots = [0]
        else:
            radical = np.sqrt(1 - norm_p ** 2)
            roots = [-radical, radical]
        self.intersection_points = [coordinates(self.p_int + self.s_int_hat * t) for t in roots]
        # self.beta = [self.gc_segment.find_beta(g_.x) for g_ in self.intersection_points]

    def arc_distance(self):
        # Calculates the distance between the great circle segments
        if not hasattr(self, 'beta'):
            self.gc_sc_intersect()
        if len(self.beta) == 1:
            d_sigma = self.gc_segment.beta_2 - self.beta[0]
        elif len(self.beta) == 2:
            d_sigma = self.beta[1] - self.beta[0]
        else:
            d_sigma = 0
        self.arc = d_sigma
        self.d = self.arc * self.gc_segment.R

def read_and_run_extreme_case():
    base_path = os.getcwd()

    #Loop over each pair of offset values, from 0.01 to 1.0e-15, stepsize 10

    start_offset = 0.01
    end_offset = 10 ** (-15)
    step_size = 10

    # Generate the offset array
    offsets = []
    value = start_offset
    while value >= end_offset:
        offsets.append(value)
        value /= step_size

    offsets = np.array(offsets)

    for i in range(len(offsets) - 1):
        start_off = offsets[i]
        end_off = offsets[i + 1]

        start_off_str = format_offset(start_off)
        end_off_str = format_offset(end_off)

        # Create filenames based on latitude ranges
        off_range = f"{start_off_str}_{end_off_str}"
        arc_file_name = f"{off_range}ArcUp_Arcs_Exponent.csv"
        arc_path = os.path.join(base_path, "generated_arcs", arc_file_name)

        # Read arcs for the current latitude range
        const_zs = []
        region_arcs, const_zs = read_great_circle_arcs_from_csv_exponent(arc_path)

        # recond a sanitized_arcs lthat will be used to store the arcs that are not empty
        sanitized_arcs = []
        sanitized_const_zs = []

        for i in range(len(region_arcs)):
            arc = region_arcs[i]
            const_z = const_zs[i]

            G_i = gc_segment(arc.get_point1(), arc.get_point2())
            S_i = sc(const_z)
            I_i = intersection(G_i, S_i)
            I_i.calculate_intersection_line()
            I_i.gc_sc_intersect()
            results = I_i.intersection_points

            # Check if the results are empty
            if len(results) > 0:
                # Check if the results not contain any NaN values or infinities
                if np.all(np.isfinite([point.x for point in results])):
                    sanitized_arcs.append(arc)
                    sanitized_const_zs.append(const_z)
            else:
                continue

        # Now rewrite the sanitized arcs to original csv file
        write_great_circle_arcs_to_csv_exponent(arc_path, sanitized_arcs, sanitized_const_zs)

        results_lat = []
        # Now loop over the sanitized arcs and calculate the intersection points
        for i in range(len(sanitized_arcs)):
            arc = sanitized_arcs[i]
            const_z = sanitized_const_zs[i]

            G_i = gc_segment(arc.get_point1(), arc.get_point2())
            S_i = sc(const_z)
            I_i = intersection(G_i, S_i)
            I_i.calculate_intersection_line()
            I_i.gc_sc_intersect()
            results = I_i.intersection_points  # there should be two intersection points

            # store the results in a list
            results_lat.append(results)

        # Now open the csv to write the results
        header = [
            "x1_significand", "x1_exponent", "y1_significand", "y1_exponent",
            "x2_significand", "x2_exponent", "y2_significand", "y2_exponent"
        ]
        file_path = os.path.join(base_path, "benchmark_results", f"{off_range}_krumm_double_ArcUp.csv")
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        with open(file_path, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)

            # Iterate through each set of results
            for result in results_lat:
                if len(result) == 2:
                    point1, point2 = result

                    # Write both points in the same row
                    writer.writerow([
                        # Point 1 (x, y)
                        DecomposedFloat(point1.x[0]).significand, DecomposedFloat(point1.x[0]).exponent,
                        DecomposedFloat(point1.x[1]).significand, DecomposedFloat(point1.x[1]).exponent,

                        # Point 2 (x, y)
                        DecomposedFloat(point2.x[0]).significand, DecomposedFloat(point2.x[0]).exponent,
                        DecomposedFloat(point2.x[1]).significand, DecomposedFloat(point2.x[1]).exponent

                    ])

                else:
                    # repeat the same point twice
                    point = result[0]
                    writer.writerow([
                        # Point 1 (x, y)
                        DecomposedFloat(point.x[0]).significand, DecomposedFloat(point.x[0]).exponent,
                        DecomposedFloat(point.x[1]).significand, DecomposedFloat(point.x[1]).exponent,

                        # Point 2 (x, y)
                        DecomposedFloat(point.x[0]).significand, DecomposedFloat(point.x[0]).exponent,
                        DecomposedFloat(point.x[1]).significand, DecomposedFloat(point.x[1]).exponent
                    ])


if __name__ == "__main__":
    read_and_run_extreme_case()
    # Define latitudes array
    latitudes = [0,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,20,30,40,50,60,70,80,89,89.9,89.99,89.990200,89.990400,89.990600,89.990800,89.999,89.9990400,89.9990800,89.9991200,89.9991600,89.9999,89.99990800,89.99991600,89.99992400,89.99993200,89.99999,90]

    #set the base path as the current directory
    base_path = os.getcwd()

    # Loop over each pair of latitude values
    for i in range(len(latitudes) - 1):
        start_lat = latitudes[i]
        end_lat = latitudes[i + 1]

        # Format latitude strings without decimal points
        start_lat_str = format_latitude(start_lat)
        end_lat_str = format_latitude(end_lat)

        # Create filenames based on latitude ranges
        lat_range = f"{start_lat_str}_{end_lat_str}"
        arc_file_name = f"{lat_range}Arcs_Exponent.csv"
        arc_path = os.path.join(base_path, "generated_arcs", arc_file_name)

        # Read arcs for the current latitude range
        const_zs = []
        region_arcs, const_zs = read_great_circle_arcs_from_csv_exponent(arc_path)

        # recond a sanitized_arcs lthat will be used to store the arcs that are not empty
        sanitized_arcs =[]
        sanitized_const_zs = []

        for i in range(len(region_arcs)):
            arc = region_arcs[i]
            const_z = const_zs[i]

            G_i = gc_segment(arc.get_point1(),arc.get_point2())
            S_i = sc(const_z)
            I_i = intersection(G_i, S_i)
            I_i.calculate_intersection_line()
            I_i.gc_sc_intersect()
            results = I_i.intersection_points

            # Check if the results are empty
            if len(results) > 0:
                # Check if the results not contain any NaN values or infinities
                if np.all(np.isfinite([point.x for point in results])):
                    sanitized_arcs.append(arc)
                    sanitized_const_zs.append(const_z)
            else:
                continue

        # Now rewrite the sanitized arcs to original csv file
        write_great_circle_arcs_to_csv_exponent(arc_path, sanitized_arcs, sanitized_const_zs)


        results_lat = []
        # Now loop over the sanitized arcs and calculate the intersection points
        for i in range(len(sanitized_arcs)):
            arc = sanitized_arcs[i]
            const_z = sanitized_const_zs[i]

            G_i = gc_segment(arc.get_point1(),arc.get_point2())
            S_i = sc(const_z)
            I_i = intersection(G_i, S_i)
            I_i.calculate_intersection_line()
            I_i.gc_sc_intersect()
            results = I_i.intersection_points # there should be two intersection points

            # store the results in a list
            results_lat.append(results)




        #Now open the csv to write the results
        header = [
            "x1_significand", "x1_exponent", "y1_significand", "y1_exponent",
            "x2_significand", "x2_exponent", "y2_significand", "y2_exponent"
        ]
        file_path = os.path.join(base_path, "benchmark_results", f"{lat_range}_krumm_double.csv")
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        with open(file_path, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)

            # Iterate through each set of results
            for result in results_lat:
                if len(result) == 2:
                    point1, point2 = result

                    # Write both points in the same row
                    writer.writerow([
                        # Point 1 (x, y)
                        DecomposedFloat(point1.x[0]).significand, DecomposedFloat(point1.x[0]).exponent,
                        DecomposedFloat(point1.x[1]).significand, DecomposedFloat(point1.x[1]).exponent,

                        # Point 2 (x, y)
                        DecomposedFloat(point2.x[0]).significand, DecomposedFloat(point2.x[0]).exponent,
                        DecomposedFloat(point2.x[1]).significand, DecomposedFloat(point2.x[1]).exponent

                    ])

                else:
                    # repeat the same point twice
                    point = result[0]
                    writer.writerow([
                        # Point 1 (x, y)
                        DecomposedFloat(point.x[0]).significand, DecomposedFloat(point.x[0]).exponent,
                        DecomposedFloat(point.x[1]).significand, DecomposedFloat(point.x[1]).exponent,

                        # Point 2 (x, y)
                        DecomposedFloat(point.x[0]).significand, DecomposedFloat(point.x[0]).exponent,
                        DecomposedFloat(point.x[1]).significand, DecomposedFloat(point.x[1]).exponent
                    ])
