
clear all

% ......................................................................
% .......... Get 3D surface points on sphere ...........................
% ......................................................................

% Author: Julia Eckert (j.eckert@imb.uq.edu.au)
% Date: 2025-5

% --- DESCRIPTION: 
% This code uses the 2D boundary of each mask per plane to identify the 
% center point of the sphere and generates the surface points and 
% corresponding normal vectors using the 'icosphere()' mesh.

% --- TO DO: 
% ImageJ/Fiji
%   Mask excl. full background
%       - create black and white mask of each image plane 
%         (background = black) ... e.g. only cells forming a lumen are
%             white, empty space in the center is black
%       - save z-stack as 'Mask_Boundary.tif'         

% Additional codes required:
%   Icosphere()     https://au.mathworks.com/matlabcentral/fileexchange/50105-icosphere

% --- INPUT:
% dimension of z-stack images, in um/pix
xyz = [0.1507813, 0.1507422, 0.4982639]; % Aggregate 

% xyz distance between planes / xyz surface points for slicing, in um 
dis_point = 6; 

% radius of the spherical system, in um
r = 40; 


% --- OUTPUT:
% save as SurfacePoints.mat
% SurfacePoints.xyz             ... xyz surface point positions on surface, in um 
% SurfacePoints.xyzNormal       ... xyz normal vectors 
% SurfacePoints.Pixel           ... xyz dimension of z-stack, in um/pix
% SurfacePoints.GridDistance    ... distance between xyz surface points, in um
% SurfacePoints.SphereRadius    ... radius of the sphere, in um
% SurfacePoints.SphereN         ... icosphere number




%..........................................................................
%................................ Main ....................................
%..........................................................................

% ... Upload images .......................................................

% z-stack image mask, background = black, aggregate = white 
image_mask = tiffreadVolume('Mask_Boundary.tif');
image_mask = double(image_mask);


% get center of mass of cell system in xy-coordinates
sm  = regionprops(bwconncomp(imfill(image_mask(:,:,round(size(image_mask,3)/2)),'holes'),8),'Area','FilledArea','FilledImage','Centroid','MajorAxisLength',...
   'MinorAxisLength', 'Perimeter','Orientation','Eccentricity','PixelIdxList','PixelList');

centroids_=[sm.Centroid];
centroids_=[centroids_(1:2:end)',centroids_(2:2:end)'];

areas_ = [sm.Area];
idx_area = find(areas_==max(areas_));

% xyz center point of cell system, in um
x_sphere_center = centroids_(idx_area,1)*xyz(1);
y_sphere_center = centroids_(idx_area,2)*xyz(2);
z_sphere_center = size(image_mask,3)/2*xyz(3);


% ... Get surface points and normals ......................................

% number of vertices based on icoshere() function: (:,1) corresponds to the
% input of icosphere() and (:,2) the number vertices as the output
Nico = [0, 12; 1, 42; 2, 162; 3, 642; 4, 2562; 5, 10242];

% surface area of the sphere with given radius
Asphere = 4*pi*(r)^2;
% approximate the number for points on sphere with distance of dis_point
Npoint = (Asphere/((3*sqrt(3)/2) * dis_point*dis_point))*3;

% get best index for icosphere mesh
[~,idx] = min(abs(Nico(:,2)-Npoint));

f = icosphere(Nico(idx,1));

% normal vectors
xyz_normal = f.Vertices;

% xyz position, in um
xyz_point = xyz_normal*r + [x_sphere_center, y_sphere_center, z_sphere_center];

SurfacePoints.xyz = xyz_point;          % xyz surface point position on surface, in um 
SurfacePoints.xyzNormal = xyz_normal;   % xyz normal vectors
SurfacePoints.Pixel = xyz;              % xyz dimension of z-stack, in um/pix
SurfacePoints.GridDistance = dis_point; % distance between xyz surface points, in um
SurfacePoints.SphereRadius = r;         % radius of the sphere, in um
SurfacePoints.SphereN = Nico(idx,1);    % icosphere number

save('SurfacePointsSphere.mat','SurfacePoints','-double')

clear all



