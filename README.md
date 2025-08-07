# GRIDGEN
A Program for Generating Unstructured Finite-Volume Grids

## Overview of GRIDGEN
GRIDGEN is a computer program for creating layered quadtree grids for use with numerical models, such as the MODFLOW-USG program for simulation of groundwater flow.

The program begins by reading a three-dimensional base grid, which can have variable row and column widths and spatially variable cell top and bottom elevations. From this base grid, GRIDGEN will continuously divide into four any cell intersecting user-provided refinement features (points, lines, and polygons) until the desired level of refinement is reached. GRIDGEN will then smooth, or balance, the grid so that no two adjacent cells, including overlying and underlying cells, differ by more than a user-specified level tolerance. Once these gridding processes are completed, GRIDGEN saves a tree structure file so that the layered quadtree grid can be quickly reconstructed as needed. Once a tree structure file has been created, GRIDGEN can then be used to:

- export the layered quadtree grid as a shapefile, 
- export grid connectivity and cell information as ASCII text files for use with MODFLOW-USG or other numerical models, and 
- intersect the grid with shapefiles of points, lines, or polygons, and save intersection output as ASCII text files and shapefiles.

The GRIDGEN program is demonstrated by creating a layered quadtree grid for the Biscayne aquifer in Miami-Dade County, Florida, using hydrologic features to control where refinement is added.

![GRIDGEN model grid from the cover of the GRIDGEN manual.](https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/side_image/public/thumbnails/image/GRIDGEN-cover-grid.jpg?itok=6YJgafTS)

## Documentation for GRIDGEN
Lien, Jyh-Ming, Liu, Guilin, and Langevin, C.D., 2015, GRIDGEN version 1.0 -- A computer program for generating unstructured finite-volume grids: U.S. Geological Survey Open-File Report 2014-1109, 39 p., http://dx.doi.org/10.3133/ofr20141109.

This report describes the theory and input instructions at the time of the initial GRIDGEN release.

## How to Cite GRIDGEN
This USGS software has two citations associated with it.

The report citation is for the original report or article documenting the underlying theory, methods, instructions, and (or) applications at the time the initial version of the software was released. This digital object identifier (DOI) is for the report.
The software release citation is for the software/code itself (now referred to by USGS as a "Software Release") and references a specific version of the code and associated release date. This DOI links to the code.
In instances where an author is citing use of this software, it would be appropriate to cite both the report documenting the code and the specific software release version that was used.

### Report Citation for GRIDGEN

Lien, Jyh-Ming, Liu, Guilin, and Langevin, C.D., 2015, GRIDGEN version 1.0 -- A computer program for generating unstructured finite-volume grids: U.S. Geological Survey Open-File Report 2014-1109, 39 p., http://dx.doi.org/10.3133/ofr20141109.

### Software/Code Citation for GRIDGEN v1.0.02

Lien, Jyh-Ming, Liu, Guilin, and Langevin, C.D., 2017, GRIDGEN version 1.0.02 -- A computer program for generating unstructured finite-volume grids: U.S. Geological Survey Software Release, 06 January 2017, http://dx.doi.org/10.5066/F79G5JXV.
