%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example code of using the DBT reconstruction functions for non-arc DBT geometry
% The projection data in this code was created by "gen_dbtproj_example.m".
%
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Load in the projection views.
%"g_noi.mat" is for SART reconstruction.
%"proj_noi.mat" is for ML reconstruction.
load example_data/g_noi.mat; %the variable is 'g_noi'.
g=g_noi;
load example_data/proj_noi.mat; %the variable is 'g_noi'.
proj=proj_noi;

%==================================
%User Defines the scanner geometry
%==================================
 dso = 60.5; %in cm: dist. from the src to the rotation ctr
 dod = 4.5;  %in cm: dist. from the rotation ctr to the detector
 dsd = dso+dod; %in cm: dist. from src to the detector
 orbit = 30; %Angular span
 na = size(g,3); %number of projection views
 ds = 0.04; %in cm: detector element pixel size; 
            %'s', the x-ray tube moving direction, positve pointing toward right. 
            % and 't', the perpendicular direction to 's' direction, positive pointing toward the nipple. 
 dt = 0.04;
 ns = size(g,1); %number of detector elements in the 's' direction
 nt = size(g,2); %...................................'t'..........
 projview_offset_s = 0; %detector center offset along the 's' direction in pixels relative to the x-ray source center
 projview_offset_t = 0; % detector center offset along the 't' direction in pixels relative to the x-ray source center
 projview_offset = [projview_offset_s projview_offset_t];


theta = [-(na-1)/2: (na-1)/2]*orbit/(na-1);  

ds = 0.04; %in cm
dt = ds;
det_depth = size(g,2)*dt;
ns = size(g,1);
nt = size(g,2);

%src_vec: a nz x 3 matrix define the x-ray source coordinates for each
%view.
%x axis: posive direction points toward right, 
%y axis: positive direction points toward the nipple
%z axis: positive direction points toward the x-ray source

%For this example, the projection date were created from a arc trajectory
%so we calculated the soruce coordinates for the projection data.
%The center of detector was treated as the origin of the 3D coordinate
%system (see readme document)
src_vec = zeros(na,3);
src_vec(:,1) = -dso*sin(theta*pi/180);
src_vec(:,2) = det_depth/2;
src_vec(:,3) = dod+dso*cos(theta*pi/180);

%=======================================
%User defines the recon volume geometry:
%=======================================
% voxel size (drx, dry, drz) in cm, 
% dimensions (nrx, nry, nrz) and 
%FOV center offsets (offset_x, offset_y, offset_z): in pixels relative to the detector center. 
%For example, if the coordinates of the FOV center is (xctr, yctr, zctr),
%then offset_x=-xctr, offset_y=-yctr and offset_z=-zctr.
nrx=504;
nry=156;
drx=0.04; 
dry=drx; 
drz=0.08; 
nrz=64;
d_objbottom_det=0; %in cm, the distance from the object bottom to the detector.
offset_x = 0;
offset_y = (det_depth/2/dry-nry/2);
zfov = (nrz-1)*drz;
offset_z = -(zfov/2+d_objbottom_det)/drz; %in pixels: The provided value here is when the object was placed right above the detector.

 %===================
 %Reconstruction
 %===================

 %Generate the system matrix
                          
 igr = image_geom('nx', nrx, 'ny',nry, 'nz', nrz, 'dx',drx, 'dz', drz,...
       'offset_x', offset_x, 'offset_y', offset_y,'offset_z', offset_z,'down', 1); 
 
 dgeo.ns = ns;
 dgeo.nt = nt;
 dgeo.ds = ds;
 dgeo.dt = dt;
     
 Gtr = Gtomo_syn_noncircular(src_vec,igr, dgeo, projview_offset);


% SART reconstruction
xbp1 = BP(Gtr, g); %initialization for SART
disp 'SART...'
[xartt1,costart] = SART_dbt(Gtr,g,xbp,2,0.5);
% ML recosntruction
disp 'ML...'
[xmlt1,costml] = ML_dbt(Gtr,proj,xbp,I0,3,2);
disp 'Recon completed';

figure('Name','dbt_recon_noncircular_example');
imagesc(xartt1(end:-1:1,:,30)), daspect([1 1 1]), colormap(gray)
title 'Slice 30 (lesion focal plane) of SART reconstruction'
colorbar

save sart_noncir.mat xartt1;
save ml_noncir.mat xmlt1; 

