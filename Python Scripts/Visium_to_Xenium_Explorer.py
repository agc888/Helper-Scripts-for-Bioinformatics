import sopa
import sopa.io
import spatialdata 
from spatialdata import SpatialData
import spatialdata_plot
import xarray as xr
#from spatialdata import SpatialImage
import anndata
import scanpy as sc
import pandas as pd
from xarray import DataArray
from spatialdata.models import Image2DModel, ShapesModel, TableModel
from spatialdata.transformations.transformations import Identity, Scale
import spatialdata_io
from shapely import geometry
import numpy as np

### Things for you to change 
path = ".../outs/"
dataset_id = "test"
explorer_path = "../Downloads/xenium_explorer_folder"


spatial_d = spatialdata_io.visium(path, dataset_id = dataset_id)


import math

# This function gets just one pair of coordinates based on the angle theta
def get_circle_coord(theta, x_center, y_center, radius):
    x = radius * math.cos(theta) + x_center
    y = radius * math.sin(theta) + y_center
    return (x,y)

# This function gets all the pairs of coordinates
def get_all_circle_coords(x_center, y_center, radius, n_points):
    thetas = [i/n_points * math.tau for i in range(n_points)]
    circle_coords = [get_circle_coord(theta, x_center, y_center, radius) for theta in thetas]
    return circle_coords

spatial_d[dataset_id]["geometry"] = spatial_d[dataset_id].apply(lambda x: geometry.Polygon(get_all_circle_coords(x["geometry"].x,
                                                    x["geometry"].y,
                                                    x["radius"],
                                                    n_points=40)), axis=1).values


sopa.io.write(
        explorer_path,
        spatial_d,
        image_key= str(dataset_id) + "_hires_image",
        ram_threshold_gb=110,
    )


#### To add metadata for plotting in xenium explorer
#### This is just an example to show how to add metadata to the explorer file
#### In this case we just asign a random number between 1-6 to each spot and call this column 'random_numbers'
#### NOTE: the metadata slot in the .obs must be astype("category") to be exported to the xenium.explorer file 


num_cells = spatial_d["table"].n_obs

# Step 2: Generate random numbers between 1 and 5 (inclusive)
random_numbers = np.random.randint(1, 6, size=num_cells)

# Step 3: Add the random numbers as a new column in the .obs slot
spatial_d["table"].obs['random_numbers'] = random_numbers
spatial_d["table"].obs['random_numbers'] = spatial_d["table"].obs["random_numbers"].astype("category")

sopa.io.write_cell_categories(explorer_path, spatial_d["table"])
