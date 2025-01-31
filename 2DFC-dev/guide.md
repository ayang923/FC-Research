# How to Use This Code

> A comprehensive guide to understand, set up, and run the code in this repository.

---

## Table of Contents

- Overview
- Setup Instructions
  - Prerequisites
  - Installation
- Code Structure
- How to Run
  - Input Details
  - Expected Outputs
- Customization
- Troubleshooting

---

## Overview

## Setup Instructions
All of this code was written and tested on MATLAB2021a but should work on other versions as well.

### Prerequisites
N/A

### Installation
N/A

## Code Structure
`FC_2D.m` contains "main" 2DFC function  
`Curve_seq_obj.m` contains class definition for "curve sequence" objects that encodes the domain of the boundary and discretization parameters  
  
`S_patch_obj.m`, `C1_patch_obj.m`, and `C2_patch_obj.m` contain class definitions for each patch type. S-type patches denote patches constructed along smooth portions of boundary (twice continuously differentiable), C1-type patches denote patches constructed along concave corners (interior angle > 180 degrees), and C2-type patches denote patches constructed along convex corners (interior angle < 180 degrees)  
`Q_patch_obj.m` contains class definition for square patches. The bulk of useable functions are defined in this class and each of the S-type, C1-type, and C2-type classes inherit from this class.  

`R_cartesian_mesh_obj.m` contains the class definition for the cartesian mesh that the FFT is computed on. The function values of each patch are interpolated on to this mesh to obtain the continuation function values on the Cartesian mesh.  

`boomerang_2D_FC.m`, `teardrop_2DFC.m`, and `polynomial_curve.m` contains examples of the 2D-FC algorithm for three different examples.

## How to Run


### Input Details

### Expected Outputs

## Customization

## Troubleshooting


