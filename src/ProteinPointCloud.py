# ProteinPointCloud.py
# Python file for storing functions related to extracting and analyzing point 
# clouds from protein files.

import os
import sys
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from random import sample
import math

from scipy.spatial.distance import directed_hausdorff

def atom_cloud(filename):
    with open(filename,"r") as outfile:
            data = outfile.readlines()
    points = []
    for line in data:
        if 'ATOM' not in line:
            continue
        try:
            point = tuple(float(x) for x in line.split()[10:13])
            if len(point) == 3:
                points.append(point)
        except:
            continue
    return np.array(points)

def skew(v):
    return np.array([[    0, -v[2],  v[1]], 
                     [ v[2],     0, -v[0]], 
                     [-v[1],  v[0],    0]])

def regularize_cloud(pcd):
    r, c = pcd.shape
    cloud = np.zeros((r,c))
    avg = np.mean(pcd, axis=0)
    max_index = -1
    max_dist = 0
    for i in range(r):
        cloud[i] = pcd[i] - avg
        dist = np.linalg.norm(cloud[i])
        if dist >= max_dist:
            max_dist = dist
            max_index = i
    cloud = cloud/max_dist
    
    align = np.array([1, 0, 0])
    v = np.cross(cloud[max_index], align)
    sin = np.linalg.norm(v)
    cos = np.dot(cloud[max_index], align)
    sv = skew(v)
    R = np.identity(3) + sv + (1 - cos)*(np.matmul(sv, sv))/(sin**2)
    for i in range(r):
        cloud[i] = np.matmul(R, cloud[i])
        
    return cloud

def hausdorff_distance(cloud_1, cloud_2):
    return max(directed_hausdorff(cloud_1, cloud_2)[0], 
               directed_hausdorff(cloud_2, cloud_1)[0])