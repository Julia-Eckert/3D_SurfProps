
clear all

% ......................................................................
% .... Get 3D surface points and normal vectors of arbitrarty shape.....
% ......................................................................

% Author: Julia Eckert (j.eckert@imb.uq.edu.au)
% Date: 2025-5

% --- DESCRIPTION: 
% This code uses the 2D boundary of each mask per plane to get the surface 
% points of the multicellular system. Local amounts of surface points are 
% used to create regression planes to get the normal vectors.

% --- TO DO: 
% ImageJ/Fiji
%   (1) Mask excl. full background
%           - create black and white mask of each image plane 
%             (background = black) ... e.g. only cells forming a lumen are
%             white, empty space in the center is black
%           - save z-stack as 'Mask_Boundary.tif'         
%   (2) Mask as (1) but with filled space
%           - expand the black and white mask of (1) and set the space enclosed 
%             by the cells to white ... e.g. the cells forming a lumen and the 
%             empty space is white too
%           - save z-stack as 'Mask_Full.tif'    

% --- INPUT:
% dimension of z-stack images, in um/pix
%xyz = [0.2071606, 0.2071606, 1];            % Zebrafish heart
xyz = [0.1507813, 0.1507422, 0.4982639];    % Aggregate

% distance between planes (i.e. xyz surface points) for slicing, in um 
dis_point = 6; 


% --- OUTPUT:
% save as SurfacePoints.mat
% SurfacePoints.xyz             ... xyz surface point positions on surface, in um 
% SurfacePoints.xyzNormal       ... xyz normal vectors 
% SurfacePoints.Pixel           ... xyz dimension of z-stack, in um/pix
% SurfacePoints.GridDistance    ... distance between xyz surface points, in um


%..........................................................................
%................................ Main ....................................
%..........................................................................

% ... Upload images .......................................................

% z-stack image mask, background = black, aggregate = white 
image_mask_region = tiffreadVolume('Mask_Boundary.tif');
image_mask = tiffreadVolume('Mask_Full.tif');



% ... Get surface points ..................................................

% ratio between x and z dimensions to get equal distribution of data points
% in all dimensions 
ratio_dimensions = xyz(3)/xyz(1);
point_select = dis_point/xyz(3);

% get xyz surface positions with specific xyz distance in um  
xyz_boundary_points_main = [];

% save 2D plane points of first and/or last plane 
additional_planes = [];
for k = 1:size(image_mask,3)
    % if mask starts/end randomly in z-stack
    if k-1 >= 1
        if k+1 <= size(image_mask,3)
            % check in which plane the mask appears first and last
            plane_previous = double(im2bw(image_mask(:,:,k-1)));
            plane_previous = sum(plane_previous(:));
            plane_next = double(im2bw(image_mask(:,:,k+1)));
            plane_next = sum(plane_next(:));
            plane_current = double(im2bw(image_mask(:,:,k)));
            plane_current = sum(plane_current(:));

            if plane_previous < 10 && plane_current > 10 && plane_current < size(image_mask,1)*size(image_mask,2)/100
                add_plane_yn = 1;
            elseif  plane_next < 10 && plane_current > 10 && plane_current < size(image_mask,1)*size(image_mask,2)/100
                add_plane_yn = 1;
            else 
                add_plane_yn = 0;
            end

            % for first and last plane having mask, create mash and add
            % surface points 
            if add_plane_yn == 1
                mesh_min = [1,1];
                mesh_max = size(image_mask,1:2);
        
                [mpx,mpy] = meshgrid(mesh_min(2):dis_point/xyz(2):mesh_max(2),...
                                mesh_min(1):dis_point/xyz(1):mesh_max(1));

                image_mask_k = double(image_mask(:,:,k));
                MM = image_mask_k*0;
                MM(round(mpx(:)),round(mpy(:))) = 1;
                MM(image_mask_k==0) = 0;
                [xx,yy] = find(MM);
                xyz_boundary_points_main = [xyz_boundary_points_main; yy*xyz(2), ...
                                    xx*xyz(1), k*xyz(3)*ones(length(xx),1)]; 

                additional_planes = [additional_planes; k]; 
            end
        end    
    end
    % if mask starts/ends in first/last image of z-stack
    if k == 1 || k == size(image_mask,3)
        plane_boundary = double(im2bw(image_mask(:,:,k)));
        plane_full = double(im2bw(image_mask_region(:,:,k)));
        if sum(plane_boundary(:)) == sum( plane_full(:))
            mesh_min = [1,1];
            mesh_max = size(image_mask,1:2);
    
            [mpx,mpy] = meshgrid(mesh_min(2):dis_point/xyz(2):mesh_max(2),...
                            mesh_min(1):dis_point/xyz(1):mesh_max(1));

            image_mask_k = double(image_mask(:,:,k));
            MM = image_mask_k*0;
            MM(round(mpx(:)),round(mpy(:))) = 1;
            MM(image_mask_k==0) = 0;
            [xx,yy] = find(MM);
            xyz_boundary_points_main = [xyz_boundary_points_main; yy*xyz(2), ...
                                xx*xyz(1), k*xyz(3)*ones(length(xx),1)]; 

            additional_planes = [additional_planes; k];            
        end
    end
end


% get xyz surface points with equal xyz distance, in um
k_range = [round(point_select)/2:round(point_select):size(image_mask,3), additional_planes']; 
for k = k_range

    im_plane = image_mask(:,:,k);
    im_plane = im2bw(im_plane);
    
    pixp = find(im_plane(1,:)==1);
    minp = min(pixp);
    maxp = max(pixp);    
    im_plane(1,minp:maxp) = 1;

    im_plane = imfill(im_plane, 'holes');
    B_boundary = bwboundaries(im_plane);
     
    % get xy boundary points of xy image
    if isempty(B_boundary) == 0
        xy_boundary_points = [];
        for kb = 1:size(B_boundary,1)
            if length(B_boundary{kb,1}(:,1)) > 30
    
                xy_boundary = B_boundary{kb,1};
                xy_boundary = [xy_boundary(:,2),xy_boundary(:,1)];
    
                xy_boundary(xy_boundary(:,2)==1,:) = [];
    
                if mod(k,2) == 1
                    shiftp = round(ratio_dimensions*point_select/2);
                else
                    shiftp = 1;
                end
            
                xy_boundary_points = [xy_boundary_points; xy_boundary(shiftp:round(ratio_dimensions*point_select):end,:).*xyz(1:2)];
            end
        end
    

        % get xyz surface points
        if length(xy_boundary_points) > 0 
            % remove points from open space 
            for xy_bp = 1:length(xy_boundary_points(:,1))
                xy_m = round(xy_boundary_points(:,1:2)./xyz(1:2));
    
                im_m = image_mask_region(:,:,k);
                seD = strel('disk',2);
                im_m = imdilate(im_m,seD);
    
                if im_m(xy_m(xy_bp,2),xy_m(xy_bp,1)) == 0
                    xy_boundary_points(xy_bp,:) = NaN;
                end
                
            end

            % remove NaN from list
            Y = ~isnan(xy_boundary_points);
            xy_boundary_points = xy_boundary_points(Y(:,1),:);
    
            % add to xyz boundary points 
            xyz_boundary_points_main = [xyz_boundary_points_main; xy_boundary_points, k*ones(length(xy_boundary_points(:,1)),1)*xyz(3)];
        
        end
    end
end

% remove double vectors 
double_vec = [];
for j = 1:length(xyz_boundary_points_main(:,1))
    idx_v = rangesearch(xyz_boundary_points_main ,xyz_boundary_points_main(j,:),(dis_point-dis_point/2));
    if length(idx_v{1,1}) > 1
        double_vec = [double_vec; idx_v{1,1}, (1:90-length(idx_v{1,1}))*NaN];
    end
end
double_vec_list = double_vec(:);
unique_list = unique(double_vec_list(double_vec_list>0));

for m = 1:length(unique_list)
    [row,col] = find(double_vec == unique_list(m));
    
    [max_n,idx_max] = max(double_vec(row));
    for n = 1:length(row)
        if n ~= idx_max
            double_vec(row(n),:) = NaN;
        end
    end
end
list_double = double_vec(:,2:end);
list_double = list_double(:);
dv = list_double(list_double>0);
xyz_boundary_points_main(dv,:) =[]; 



% ... Get normal vectors ..................................................

% get normal vector for all surface points via regression planes
xyz_point = [];
xyz_normal = [];
for k_point = 1:length(xyz_boundary_points_main(:,1))    

    disp([num2str(k_point), ' of ', num2str(length(xyz_boundary_points_main(:,1)))]);

    % xyz_pos of data points within a certain radius 
    Idx = rangesearch(xyz_boundary_points_main(:,1:3),xyz_boundary_points_main(k_point,1:3),dis_point*2.0);

    % position of plane points for regression
    xyz_pos = xyz_boundary_points_main(Idx{1,1},:);

    center_pos = nanmean(xyz_pos,1);
    Y = xyz_pos - center_pos;
    [~,~,V] = svd(Y,0);
     
    % normal vector of regression plane
    nP = V(:,3)';

    nP = nP*2;
    % test direction of normal vector (towards mask or away) and correct for direction 
    if round((xyz_boundary_points_main(k_point,3))/xyz(3)+nP(3)) < 1
        n_direction = 1;
    elseif round((xyz_boundary_points_main(k_point,3))/xyz(3)+nP(3)) >= size(image_mask,3)  
        n_direction = 1;
    else
        test_mask = double(image_mask_region(:,:,round((xyz_boundary_points_main(k_point,3))/xyz(3)+nP(3))));
        test_mask = im2bw(test_mask);
        test_y = round((xyz_boundary_points_main(k_point,1))/xyz(1)+nP(1));
        test_x = round((xyz_boundary_points_main(k_point,2))/xyz(2)+nP(2));
    
        if test_x <= size(image_mask,1) && test_y <= size(image_mask,2)
            if test_x > 0 && test_y > 0
                test_pix = test_mask(test_x,test_y);
        
                if test_pix == 1
                    n_direction = -1;
                else
                    n_direction = 1;
                end
            end
        else
            n_direction = 1;
        end
    end
    nP = nP/2;

    xyz_point = [xyz_point; xyz_boundary_points_main(k_point,:)];
    xyz_normal = [xyz_normal; nP*sign(n_direction)];
end

SurfacePoints.xyz = xyz_point;              % xyz position of surface points, in um
SurfacePoints.xyzNormal = xyz_normal;       % xyz normal vectors
SurfacePoints.Pixel = xyz;                  % dimension of z-stack, um/pix
SurfacePoints.GridDistance = dis_point;     % distance between xyz surface points, in um
save('SurfacePoints.mat','SurfacePoints','-double')

clear all
