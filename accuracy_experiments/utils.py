import math
import sys
import logging
from decimal import Decimal, getcontext
import os

# Setup logging
logging.basicConfig(level=logging.ERROR, stream=sys.stderr)
logger = logging.getLogger(__name__)

# Precision traits for different types
class PrecisionTraits:
    precision_map = {
        float: 24,   # 24 bits for float mantissa (single precision)
        Decimal: 53  # 53 bits for double precision, similar to C++ double
    }

    @classmethod
    def get_precision(cls, type_):
        return cls.precision_map.get(type_, None)

# Function to log errors
def log_error(method_name, message, value, significand, exponent):
    logger.error(f"Error in method: {method_name}, Value: {value}, Message: {message}, "
                 f"Significand: {significand}, Exponent: {exponent}")

class DecomposedFloat:
    def __init__(self, value=None, significand=0, exponent=0):
        self.significand = significand
        self.exponent = exponent
        if value is not None:
            self.decompose(value)

    def decompose(self, v):
        if v == 0:
            self.significand = 0
            self.exponent = 0
        else:
            # Determine the precision based on the type of the input
            precision = 53
            if precision is None:
                raise TypeError("Unsupported type for decomposition")

            # Use `math.frexp` to get the exponent and the mantissa
            mantissa, exp = math.frexp(v)
            self.exponent = exp - precision  # Adjust exponent based on precision
            value_exp = math.ldexp(mantissa, precision)
            self.significand = int(value_exp)

            # Check if the integer conversion was exact
            if float(self.significand) != value_exp:
                log_error("DecomposedFloat", "Integer part was not integer!", v, self.significand, self.exponent)
                raise ValueError("Integer part was not integer!")

            # Verify reconstruction
            if self.to_float(type(v)) != v:
                log_error("DecomposedFloat", "Could not reconstruct the original value!", v, self.significand, self.exponent)
                raise ValueError("Could not reconstruct the original value!")

    def __str__(self):
        return f"{self.significand} * 2**{self.exponent}"

    def to_float(self, type_):
        return math.ldexp(type_(self.significand), self.exponent)

    def to_csv(self):
        return f"{self.significand},{self.exponent}"

import csv

class Arc_T:
    def __init__(self):
        # Initialize a 3x2 matrix (assuming `arc(row, col)` indexing in C++ meant a 3x2 structure)
        self.data = [[0.0 for _ in range(2)] for _ in range(3)]

    def __setitem__(self, index, value):
        row, col = index
        self.data[row][col] = value

    def __getitem__(self, index):
        row, col = index
        return self.data[row][col]

    def get_point1(self):
        """Retrieve the first column (x, y, z) as a list."""
        return [self.data[row][0] for row in range(3)]

    def get_point2(self):
        """Retrieve the second column (x, y, z) as a list."""
        return [self.data[row][1] for row in range(3)]

# Helper function to reconstruct a float from significand and exponent
def reconstruct_double(significand, exponent):
    return math.ldexp(float(significand), exponent)

def read_great_circle_arcs_from_csv_exponent(filename):
    arcs = []
    const_zs = []

    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header line

        for row in reader:
            significands = []
            exponents = []

            # Parse each row alternating between significand and exponent
            toggle = True
            for value in row:
                if toggle:
                    significands.append(int(value))  # Read significand
                else:
                    exponents.append(int(value))  # Read exponent
                toggle = not toggle  # Alternate between significand and exponent

            # Reconstruct Arc_T with the parsed significands and exponents
            arc = Arc_T()
            for i in range(6):
                row_idx = i % 3
                col_idx = i // 3
                arc[row_idx, col_idx] = reconstruct_double(significands[i], exponents[i])
            arcs.append(arc)

            # Reconstruct constZ and add to const_zs
            const_z = reconstruct_double(significands[6], exponents[6])
            const_zs.append(const_z)

    return arcs, const_zs


def format_latitude(latitude):
    # Format the latitude with six decimal places
    lat_str = f"{latitude:.6f}"

    # Replace the decimal point with an underscore
    lat_str = lat_str.replace('.', '_')

    return lat_str

def format_offset(offset):
    # Format the latitude with six decimal places
    lat_str = f"{offset:.16f}"

    # Replace the decimal point with an underscore
    lat_str = lat_str.replace('.', '_')

    return lat_str


def write_great_circle_arcs_to_csv_exponent(filename, arcs, const_zs):
    """
    Function to generate and write great circle arcs to a CSV file in an exponential format.

    Parameters:
    - filename: The output CSV file name.
    - arcs: List of Arc_T objects representing great circle arcs.
    - const_zs: List of double values representing constZ values.
    """
    # Check if arcs list is empty
    if not arcs:
        print("No arcs to write to CSV.")
        return

    # Check if the file exists and delete it
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except OSError as e:
            print(f"Error deleting file: {e}")
            return  # Return on error

    # Open file in write mode
    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write CSV headers
        writer.writerow([
            "pointA_x_significand", "pointA_x_exponent", "pointA_y_significand", "pointA_y_exponent",
            "pointA_z_significand", "pointA_z_exponent", "pointB_x_significand", "pointB_x_exponent",
            "pointB_y_significand", "pointB_y_exponent", "pointB_z_significand", "pointB_z_exponent",
            "constZ_significand", "constZ_exponent"
        ])

        # Process each arc and write decomposed values
        for i, arc in enumerate(arcs):
            # Calculate constZ from const_zs list
            const_z = const_zs[i]
            df_const_z = DecomposedFloat(const_z)

            # Prepare row data for point A and point B
            row_data = []
            for col in range(2):  # 2 columns for Point A and Point B
                for row in range(3):  # x, y, z coordinates
                    df = DecomposedFloat(arc[row, col])
                    row_data.extend([df.significand, df.exponent])

            # Add constZ decomposed values to row_data
            row_data.extend([df_const_z.significand, df_const_z.exponent])

            # Write row data to CSV file
            writer.writerow(row_data)