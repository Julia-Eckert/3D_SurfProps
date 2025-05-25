
clear all

% ......................................................................
% ........................ Figures and Plots ...........................
% ......................................................................

% Author: Julia Eckert (j.eckert@imb.uq.edu.au)
% Date: 2025-5

% --- DESCRIPTION: 
% This code provides several options to visualize the nematic director
% field. Each section is seperated by '%%', starts with a small 
% description at the beginning, and is executed by Run Section. 
% These parts can be used as a template. 


%..........................................................................
%................................ Main ....................................
%..........................................................................


% To execute the section of interest, click on the section and
% press 'Run Section'


%% .........For 3D image reconstruction of the z-stack.....................

% This section reconstructs the 3D image of the multicellular system. The
% processing time depends on the size of z-stack and is specified in the
% Command Window.

disp(['Start plotting']);

load('Analysis_Surface_Nematic.mat')
load('SurfacePoints.mat', 'SurfacePoints')

image_mask = tiffreadVolume('Mask_Boundary.tif');
image_boundary = tiffreadVolume('Orientation_Ch.tif');

% change to grayscale image
tempb = [];
for km = 1:size(image_boundary,3)
    tempb(:,:,km) = imadjust(image_boundary(:,:,km));
end
image_boundary = tempb;

% xyz dimension of the z-stack, in um/pix
xyz = SurfacePoints.Pixel;

% interpolate and remove background
xlin = 0:xyz(2):size(image_boundary,2)*xyz(2)-xyz(2);
ylin = 0:xyz(1):size(image_boundary,1)*xyz(1)-xyz(1);
zlin = 0:xyz(3):size(image_boundary,3)*xyz(3)-xyz(3);
[xdata, ydata, zdata] = meshgrid(xlin,ylin,zlin);

% changes the size of the mask slightly
image_boundary = double(image_boundary);
se = strel('disk',max(Properties.range_tangent_planes));
image_mask = imerode(image_mask,se);
%image_mask = imdilate(image_mask,se);
image_boundary(image_mask==0)=NaN;

% new image of system
vdata = double(image_boundary);
xmin = min(xdata(:));
xmax = max(xdata(:));
ymax = max(ydata(:));
ymin = min(ydata(:));
zmax = max(zdata(:));
zmin = min(zdata(:));

% plot
figure
delta_plane = 0:xyz(3):zmax;
for kk = 1:length(delta_plane)
    ks = delta_plane(kk);
    disp([num2str(kk), ' of ', num2str(length(delta_plane))]);
    %disp([num2str(ks+1), ' of ', num2str(round(zmax/xyz(3))+1)]);
    hslice = surf(linspace(xmin,xmax,1000), linspace(ymin,ymax,1000), ones(1000)*ks);
    rotate(hslice,[-1,0,0],0) % rotation axis and angle
    xda = get(hslice,'XData')-Properties.x_system_center;
    yda = get(hslice,'YData')-Properties.y_system_center;
    zda = get(hslice,'ZData')-Properties.z_system_center;
    delete(hslice)
    %show slice of interest
    hold on
    h = slice(xdata-Properties.x_system_center, ydata-Properties.y_system_center, zdata-Properties.z_system_center, vdata, xda, yda, zda);
    set(h, 'EdgeColor','none', 'FaceColor','interp', 'DiffuseStrength', 0.8)
    alpha(.9)
end
colormap('gray');
axis equal

disp(['End plotting']);



%% ..........For plotting nematic director field and defects................

% This section plots the nematic director field, where the directors are 
% black. The field is overlayed with nematic defects.


disp(['Start plotting']);

load('Analysis_Coarse_Grained_Nematic.mat')
load('SurfacePoints.mat')
load('Analysis_Surface_Nematic.mat')

xyz_pos = SurfacePoints.xyz;
Nem_director = CoarseGrainedNematic.Nematic_Director;
Q = CoarseGrainedNematic.Local_Order;
scale_factor = 3;

x_system_center = Properties.x_system_center;
y_system_center = Properties.y_system_center;
z_system_center = Properties.z_system_center;

figure
hold on
% plot director field, color = black
for k = 1:length(xyz_pos(:,1))
    U = Nem_director(k,:);
    hDefl = quiver3(xyz_pos(k,1)-x_system_center,xyz_pos(k,2)-y_system_center,xyz_pos(k,3)-z_system_center,-scale_factor*U(1),-scale_factor*U(2),-scale_factor*U(3));
    set(hDefl,...
    'Color','k',...
    'LineWidth',3,...
    'MaxHeadSize',0,...
    'AutoScale','off');
    hDefl = quiver3(xyz_pos(k,1)-x_system_center,xyz_pos(k,2)-y_system_center,xyz_pos(k,3)-z_system_center,scale_factor*U(1),scale_factor*U(2),scale_factor*U(3));
    set(hDefl,...
    'Color','k',...
    'LineWidth',3,...
    'MaxHeadSize',0,...
    'AutoScale','off');
end

% plot defects
hold on
load('Analysis_Defects.mat')
DC_all = Defects.Defect_Charge;
xyz_pos_defect = Defects.xyz_Defect_Charge;
for kk = 1:length(DC_all)
    hold on
    XY = xyz_pos_defect(kk,:);
    DC = DC_all(kk);
    if DC < -0.8 && DC > -1.2 %-1
        plot3(XY(1)-x_system_center, XY(2)-y_system_center,XY(3)-z_system_center,'o', 'color',[0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 2)
    end
    if DC < -0.3 && DC > -0.7 %-1/2
        plot3(XY(1)-x_system_center, XY(2)-y_system_center,XY(3)-z_system_center, 'o','color',[0.3010 0.7450 0.9330], 'MarkerFaceColor', 'cyan','MarkerSize', 15, 'LineWidth', 2)
    end
    if DC < 0.7 && DC > 0.3 %+1/2
        plot3(XY(1)-x_system_center, XY(2)-y_system_center,XY(3)-z_system_center,'o', 'color',[0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerSize', 15, 'LineWidth', 2)
    end
    if DC < 1.2 && DC > 0.8 %+1
        plot3(XY(1)-x_system_center, XY(2)-y_system_center,XY(3)-z_system_center,'o', 'color',[0.6350 0.0780 0.1840], 'MarkerSize', 15, 'LineWidth', 2)
    end
end

disp(['End plotting']);


%% ..For plotting nematic director field color-coded with nematic order................

% This section plots the nematic director field, where the directors are 
% color-coded in respect to the nematic order: red = order of 1, blue = 0.


disp(['Start plotting']);

load('Analysis_Coarse_Grained_Nematic.mat')
load('SurfacePoints.mat')
load('Analysis_Surface_Nematic.mat')

xyz_pos = SurfacePoints.xyz;
Nem_director = CoarseGrainedNematic.Nematic_Director;
Q = CoarseGrainedNematic.Local_Order;
scale_factor = 3;

x_system_center = Properties.x_system_center;
y_system_center = Properties.y_system_center;
z_system_center = Properties.z_system_center;

figure
hold on
% plot director field, red: order = 1, blue: order = 0
for k = 1:length(xyz_pos(:,1))
    U = Nem_director(k,:);
    cmap = colormap('jet');
    m = 0:1/(length(cmap(:,1))-1):1;
    idx = rangesearch(m',Q(k,1),1/(length(cmap(:,1))-1));
    if length(idx{1,1}) > 0 
        rm = cmap(idx{1,1},:);
        rm = rm(1,:);
        hDefl = quiver3(xyz_pos(k,1)-x_system_center,xyz_pos(k,2)-y_system_center,xyz_pos(k,3)-z_system_center,-scale_factor*U(1),-scale_factor*U(2),-scale_factor*U(3));
        set(hDefl,...
        'Color',rm,...
        'LineWidth',3,...
        'MaxHeadSize',0,...
        'AutoScale','off');
        hDefl = quiver3(xyz_pos(k,1)-x_system_center,xyz_pos(k,2)-y_system_center,xyz_pos(k,3)-z_system_center,scale_factor*U(1),scale_factor*U(2),scale_factor*U(3));
        set(hDefl,...
        'Color',rm,...
        'LineWidth',3,...
        'MaxHeadSize',0,...
        'AutoScale','off');
    end
end

disp(['End plotting']);