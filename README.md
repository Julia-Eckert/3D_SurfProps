## 1 Cite as
Eckert et al. (2025) Capturing Nematic Order on Tissue Surfaces of Arbitrary Geometry. bioRxiv. [https://doi.org/10.1038/s41567-023-02179-0](https://www.biorxiv.org/content/10.1101/2025.01.28.635015v1)

Author: Julia Eckert (j.eckert@imb.uq.edu.au)

Date: 2025-5

## 2 Description 

Here, we provide MATLAB scripts of our analysis pipeline that properly captures the nematic order of tissue surfaces of arbitrary geometry. The analysis is based on 2D images of a single z-stack acquired with a conventional widely available microscope. It does so by slicing the tissue surface with translated tangent planes, i.e. shallow slices, and reconstructing the order parameter field within these slices, resulting in directors that are truly tangent to the surface and do not lose information due to geometry. 

![alt text](https://github.com/Julia-Eckert/3D_SurfProps/blob/main/Examples/Figure.png)

A: Multicellular MCF10A aggregate. Dark areas indicate the actin filaments at the cell boundaries, while lighter areas of the cells are used to perform the orientation analysis using OrientationJ. The image is superimposed with the nematic director field (white lines) and nematic topological defects (orange: +1/2, cyan: -1/2). 

B: Outer compact layer of the ventricle of the zebrafish heart, 120 hpf. Dark areas indicate the fluorescence signal of the membrane (myl7:BFP-CAAX) at the cell
boundaries. The image is superimposed with the nematic directors, which are color-coded with the local nematic order parameter (red: high alignment of cells and order of 1, blue: misalignment and order of 0).


## 3 Required software

- [ImageJ/Fiji](https://imagej.net/software/fiji/)
- [OrientationJ](https://bigwww.epfl.ch/demo/orientation/) plugin in ImageJ/Fiji 
- [mij.jar](https://github.com/dasv74/mij) and [ij.jar](https://github.com/dasv74/mij) for running ImageJ/Fiji within MATLAB



## 4 Required z-stacks and image pre-processing 

To perform the orientation analysis, the scripts below require three z-stacks as input files. See examples in the *Examples* folder.  

**Orientation_Ch.tif**     
This z-stack contains the structure of interest (e.g. F-actin, membrane marker, etc.) to be used for the orientation analysis. The aim of this image pre-processing is to obtain a uniform intensity signal in each image and across the entire z-stack. Image pre-processing can be performed in Fiji/ImageJ with the following steps: (1) subtracting the background, (2) smoothening the intensity signal with the Gaussian Blur filter, (3) uniforming the intensity of each image with the Normalized Local Contrast plugin, and (4) uniforming the intensity across the z-stack with the [Stack Contrast Adjustment](https://imagej.net/ij/plugins/stack-contrast/index.htm) plugin. If necessary, the z-stack must be inverted before saving with the structure for orientation analysis in bright pixels (e.g. membrane marker as dark pixels to analyze the orientation of the cells in bright pixels). 

**Mask_Boundary.tif**     
This z-stack contains the black-and-white mask (cell system: white, background: black) of each image. Threshold methods are provided in Fiji/ImageJ. 

**Mask_Full.tif**     
This z-stack contains the black-and-white mask (cell system: white, background: black) of each image as described above, but possible empty spaces are also set to white (e.g. the cells that form a lumen and the empty space are white). Both masks can be identical and are required to determine surface points, normal vectors, and the tangent planes for the analysis. 



## 5 Scripts

Each of the following scripts contains a detailed description at the beginning of the script. This includes the input parameters required and format of the images / z-stack. 

---

### 5.1 Surface points and normal vectors

We provide two scripts to determine the xyz surface point positions with their corresponding 3D normal vectors. Other software can be used to obtain these properties. 

```MATLAB
Surface_Points.m
```
... This script is used for multicellular systems with arbitrary geometry. It creates 2D boundaries of each mask per image of the z-stack to obtain xy boundary points, which are then the corresponding xyz surface points of the multicellular system. The **distance between the surface points** is selected by the user and influences two factors. (1) Local sets of surface points are used to create regression planes to determine the normal vectors. By default, the set of surface points is set to 2 times the distance between the surface points. If the distance is set too large, this may affect the orientation of the normal vector to the surface, depending on the curvature of the system. (2) The distance between the surface points is equal to the distance between the nematic directors. To account for local variations in cell/tissue orientation and local changes in curvature, we recommend not increasing the distance beyond the cell diameter.   

**--Input:**
- **Mask_Boundary.tif** ... z-stack image mask of cell system, background = black, cell system = white 
- **Mask_Full.tif**     ... z-stack image mask of cell system with filled space, background = black, cell system = white 
- **xyz(1:3)**          ... list providing the dimension of the z-stack, in µm/pixel 
- **dis_point**         ... value that defines the distance between planes (i.e. xyz surface points) for slicing, in µm 

**--Output:**

*SurfacePoints.mat*
- SurfacePoints.xyz(k,1:3)        ... xyz position of *k* surface points on the surface, in µm
- SurfacePoints.xyzNormal(k,1:3)  ... xyz normal vectors on surface 
- SurfacePoints.Pixel(1:3)        ... xyz dimension of the z-stack (=Input), in µm/pixel
- SurfacePoints.GridDistance      ... distance between xyz surface points (=Input), in µm

...............

```MATLAB
Surface_Points_Sphere.m
```
... This script can be used for systems with a spherical shape. It determines the xyz surface points and normal vectors based on a homogeneous triangular mesh using the *[icosphere()](https://au.mathworks.com/matlabcentral/fileexchange/50105-icosphere)* function.

**--Input:**
- **Mask_Boundary.tif** ... z-stack image mask of cell system, background = black, cell system = white 
- **xyz(1:3)**          ... list providing the dimension of the z-stack, in µm/pix 
- **dis_point**         ... value that defines the xyz distance between planes (i.e. xyz surface points) for slicing, in µm 
- **r**                 ... value that defines the radius of the spherical multicellular system used for the **icosphere()** function, in µm

**--Output:**

*SurfacePointsSphere.mat*
- SurfacePoints.xyz(k,1:3)        ... xyz position of *k* surface points on the surface, in µm
- SurfacePoints.xyzNormal(k,1:3)  ... xyz normal vectors on surface 
- SurfacePoints.Pixel(1:3)        ... xyz dimension of the z-stack (=Input), in µm/pixel
- SurfacePoints.GridDistance      ... distance between xyz surface points (=Input), in µm
- SurfacePoints.SphereRadius      ... radius of the sphere (=Input), in µm
- SurfacePoints.SphereN           ... value of N subdivisions as defined by the **icosphere()** function

---

### 5.2 Tangent planes and analysis of shallow slices

The following scripts generate tangent planes located at the xyz surface points and translated along the normal vectors in the direction of the bulk, creating a shallow slice onto which the signal of interest is projected. Only local information of these slices is taken and projected back onto the surface of the 3D system. The **distance of the slices to the surface** and whether the signal is **max or mean projected over several slices** are selected by the user. 


```MATLAB
Slicing_Nematic.m
```
... This script determines the nematic orientation field using the plugin [OrientationJ](https://bigwww.epfl.ch/demo/orientation/) in ImageJ/Fiji. It requires a connection between ImageJ and MATLAB via [MIJ](https://imagej.net/plugins/miji). This connection is ensured by the Miji.m script, which is automatically provided by the OrientationJ plugin. To successfully run the script, the user is asked to install the required software as listed in **3 Required software**. In order for MATLAB to access these files, the paths in the script need to be changed. If the connection between MATLAB and ImageJ/Fiji is successful, the ImageJ/Fiji window will open when the script is executed and show the orientation analysis. If ImageJ/Fiji opens, but an error message appears, the connection to Miji.m was not successful and needs to be checked again. In addition, the command window shows the progress of executing the script by displaying the successfully analyzed surface points. For example: 1 of 80. The analysis of one surface point / tangent plane can take approximately 45s. 


**--Input:**
- **Orientation_Ch.tif**    ... z-stack of signal for orientation analysis, structure of interest = bright pixels
- **Mask_Boundary.tif**     ... z-stack image mask of cell system, background = black, cell system = white 
- **SurfacePoints.mat**     ... output file generated by Surface_Points.m, see above 5.1
- **no_tang_plane**         ... value or list of slice distances per surface point to the surface, e.g. [1:1:3], in µm
- **cg_range**              ... value for coarse-graining range: radius = cg_range * W with W = dis_point/2 (see article, Eq.(4))

**--Output:**

*Analysis_Surface_Nematic.mat*

Plane analysis:
- Nematic(k).xyz_director_pos(1:3)         ... xyz position of the *k* surface point, in µm
- Nematic(k).xyz_normal_vector(1:3)        ... xyz normal vector
- Nematic(k).theta                         ... theta angle of plane, in degree
- Nematic(k).phi                           ... phi angle of plane, in degree
- Nematic(k).Nematic_director_comp         ... 2D complex nematic order 
- Nematic(k).xy_Nematic_director           ... xy nematic director from the tangent plane analysis

Properties:
- Properties.x_system_center               ... x cell system center point, in µm
- Properties.y_system_center               ... y cell system center point, in µm
- Properties.z_system_center               ... z cell system center point, in µm
- Properties.grid_dis_3D                   ... distance between xyz surface points, in µm
- Properties.range_tangent_planes          ... distance of slices to surface (=Input), e.g. [1:1:3], in µm
- Properties.coarse_graining_range         ... value of coarse-graining range (=Input)

...............

```MATLAB
Slicing_Intensity.m
```
... This script determines the mean and max intensity near the surface point. 

**--Input:**
- **Marker_X.tif**          ... non pre-processed z-stack of the fluorescent marker
- **Mask_Boundary.tif**     ... z-stack image mask of cell system, background = black, cell system = white 
- **SurfacePoints.mat**     ... output file generated by Surface_Points.m, see above 5.1
- **no_tang_plane**         ... value or list of slice distances per surface point to the surface, e.g. [1:1:3], in µm
- **r_range**               ... value for coarse-graining range: radius = r_range * dis_point 

**--Output:**

*Analysis_Surface_Marker_X.mat*

Background signal of the *k*-th tangent slice:
- X_Intensity.Bkg_max(k,1:3)           ... (k,1) median, (k,2) mean, (k,3) s.d. of max-projected signal  
- X_Intensity.Bkg_sum(k,1:3)           ... (k,1) median, (k,2) mean, (k,3) s.d. of sum-projected signal 
- X_Intensity.Bkg_mean(k,1:3)          ... (k,1) median, (k,2) mean, (k,3) s.d. of mean-projected signal

Fluorescence signal within the coarse-graining radius of the *k*-th surface point: 
- X_Intensity.max_projection_mean(k,1:3)    ... (k,1) median, (k,2) mean, (k,3) s.d. of max-projected signal 
- X_Intensity.sum_projection_mean(k,1:3)    ... (k,1) median, (k,2) mean, (k,3) s.d. of sum-projected signal 
- X_Intensity.mean_projection_mean(k,1:3)   ... (k,1) median, (k,2) mean, (k,3) s.d. of mean-projected signal 
- X_Intensity.radius                        ... radius of region of interest, in µm

Properties:
- Properties.x_system_center               ... x cell system center point, in µm
- Properties.y_system_center               ... y cell system center point, in µm
- Properties.z_system_center               ... z cell system center point, in µm
- Properties.grid_dis_3D                   ... distance between xyz surface points, in µm
- Properties.range_tangent_planes          ... distance of slices to surface (=Input), e.g. [1:1:3], in µm
- Properties.coarse_graining_range         ... value of coarse-graining range (=Input)

---

### 5.3 Nematic director field and topological defects

```MATLAB
Nematic_Defect.m
```
... This script coarse-grains the nematic orientation field on 3D surfaces, calculates the local nematic order parameter, and detects nematic topological defects. The **coarse-graining radius** is determined by the user. cg_range = 0 is the non-coarse-grained field. The larger the radius, the smoother the orientation field becomes, and nearby defect pairs cancel each other out. The script presented in **5.4 Plots and figures** provides an option to visualize the director field.

**--Input:**
- **Analysis_Surface_Nematic.mat**     ... output file generated by Slicing_Nematic.m, see above 5.2
- **cg_range**                         ... value for coarse-graining range, r = cg_range * Properties.grid_dis_3D 
- **analysis_defects**                 ... = 1: for defect analysis, = 0: for no defect analysis 
- **analysis_sphere**                  ... = 1: for a sphere, = 0 others
   
**--Output:**

*Analysis_Coarse_Grained_Nematic.mat*

- CoarseGrainedNematic.Nematic_Director(k,1:3)       ... xyz coarse-grained nematic director at the *k*-th surface point
- CoarseGrainedNematic.Local_Order(k,1)              ... local nematic order parameter of coarse-grained director field
- CoarseGrainedNematic.Coarse_Graining_radius        ... coarse-grained radius, in µm
- CoarseGrainedNematic.Coarse_Graining_number_points ... local number of surface points for coarse-graining

*Analysis_Defects.mat*

- Defects.xyz_Defect_Charge(m,1:3)                 ... xyz position of the *m*-th defect, in µm
- Defects.Defect_Charge(m)                         ... defect charge    

---

### 5.4 Plots and figures

The output files in 5.1 to 5.3 contain all information to visualize the nematic director field with topological defects, the local nematic order parameter, and the intensity profile. Here, we provide simple scripts to visualize some properties. Each section separated by %% can be run individually.

```MATLAB
Plots.m
```


## 6 Examples

The *Examples* folder contains two analyzed multicellular systems (see **2 Description**). Each of the folders contains the MATLAB output files and a folder named  *Single_image*, which includes an example of a single plane of the z-stack of each input file to show what the black-and-white masks and the channel for the orientation analysis look like. The full z-stack can be found on Zenodo.

### 6.1 Multicellular aggregate

Input for *Surface_Points.m*:
```MATLAB
xyz = [0.1507813, 0.1507422, 0.4982639];
dis_point = 6;
```

Input for *Slicing_Nematic.m*:
```MATLAB
no_tang_plane = 5;
cg_range = 2.0;
```

Input for *Nematic_Defect.m*:
```MATLAB
cg_range = 1.7;
```

### 6.2 Zebrafish heart

Input for *Surface_Points.m*:
```MATLAB
xyz = [0.2071606, 0.2071606, 1];   
dis_point = 10;
```

Input for *Slicing_Nematic.m*:
```MATLAB
no_tang_plane = [0:1:3];
cg_range = 2.0;
```

Input for *Analysis_Surface_Marker_X.m*:
```MATLAB
no_tang_plane = [0:1:3];
r_range = 0.5;    
```

Input for *Nematic_Defect.m*:
```MATLAB
cg_range = 1.7;
```
