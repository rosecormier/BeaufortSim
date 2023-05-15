import sys
import argparse
import matplotlib.pyplot as plt

from ecco_visualization import *

###############

#Parse command-line input

parser = argparse.ArgumentParser(description="Plot pressure and velocity fields in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, default=[65.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, default=[90.0, 180.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 4])
parser.add_argument("--datdir", type=str, help="Directory (relative to home) to store ECCO data", default="/Downloads/")
parser.add_argument("--outdir", type=str, help="Output directory (relative to here)", default="/visualization/")

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()

config = vars(args)
print(config)