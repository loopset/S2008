import numpy as np
from numpy.typing import NDArray
import awkward as ak
import pandas as pd
import numpy as np
from tensorflow import keras
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
import ROOT

# Read labels
dflabel = pd.read_csv("./events_f0.csv")
last_valid = dflabel["type"].last_valid_index()

# Encode labels
le = LabelEncoder()
y = le.fit_transform(dflabel.iloc[:last_valid]["type"])


ROOT.EnableImplicitMT()

rdf = ROOT.RDataFrame("SimpleTree", ["./Outputs/simple_tree.root"])  # type: ignore

arr = ak.from_rdataframe(rdf, columns=("X", "Y", "Z", "Ef0", "fRun", "fEntry"))


def prepare_spatial(events: ak.Array, rebin_factor=2, max_coord=128):
    """
    Project 3D hits to XY plane and rebin into binary occupancy grid.
    Drops Z, deduplicates (X,Y) pairs, rebins, returns (N_events, grid_size, grid_size, 1).

    Args:
        events: awkward array (N_events, variable_hits, 3)
        rebin_factor: bin aggregation factor
        max_coord: maximum coordinate value in input (e.g., 128 for [0, 127])
    """
    grid_size = max_coord // rebin_factor
    grids = np.zeros((len(events), grid_size, grid_size, 1), dtype=np.float32)

    for i, event in enumerate(events):
        # X
        x = ak.to_numpy(event["X"]).astype(int)
        # Y
        y = ak.to_numpy(event["Y"]).astype(int)
        xy = np.column_stack([x, y])
        # Rebin
        xy_rebinned = (xy // rebin_factor).astype(int)
        # Find unique
        unique_xy = np.unique(xy_rebinned, axis=0)
        unique_xy = unique_xy[
            (unique_xy >= 0).all(axis=1) & (unique_xy < grid_size).all(axis=1)
        ]
        grids[i, unique_xy[:, 0], unique_xy[:, 1], 0] = 1

    return grids


X_space = prepare_spatial(arr, rebin_factor=4)
X_scalar = ak.to_numpy(arr["Ef0"])

# Labels
