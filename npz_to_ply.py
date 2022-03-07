import numpy as np
import open3d as o3d

file="output_data/ETH3D/pipes/raw/lidar.npz" # "lidar.npz" or "mvs.npz"
data=np.load(file)
points = data["points"]
sensor_position = data["sensor_position"]

pcd = o3d.geometry.PointCloud()   # create an empty point cloud object
pcd.points = o3d.utility.Vector3dVector(data["points"].astype(np.float16))   # fill with points
pcd.normals = o3d.utility.Vector3dVector(data["sensor_position"])   # fill normal field with sensor position
o3d.io.write_point_cloud(file[:-4]+".ply", pcd)    # export point cloud
