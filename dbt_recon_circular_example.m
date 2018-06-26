%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example code of using the DBT reconstruction functions for arc DBT geometry.
% The projection data in this code was created by "gen_dbtproj_example.m".
%
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%load in the projection views.
%"g_noi.mat" is for FBP and SART reconstruction.
%"proj_noi.mat" is for ML reconstruction.
load example_data/proj_noi.mat; %the variable is 'proj_noi'.
load example_data/g_noi.mat; %the variable is 'g_noi'.
g=g_noi;
proj=proj_noi;

% #### If Breast is in the right side ####
% g = flip(rot90(g,2),3);
% proj = flip(rot90(proj,2),3);


%==================================
%User Defines the scanner geometry
%==================================
 
 dso = 60.5; %in cm: dist. from the source to the rotation center
 dod = 4.5;  %in cm: dist. from the rotation center to the detector
 dsd = dso+dod; %in cm: dist. from source to the detector
 
 orbit =30;  %in degree: angular span
 na = size(g,3); %number of projection views
 ds = 0.04; %in cm: detector element pixel size in the 's' direction; 
 dt = 0.04; %in cm: detector element pixel size in the 's' direction; 
            %'s', the x-ray tube moving direction, positive pointing toward right.           
            %'t', the perpendicular direction to 's' direction, positive pointing toward the nipple.            
 
 ns = size(g,1); %number of detector elements in the 's' direction
 nt = size(g,2); %number of detector elements in the 's' direction
 offset_s = 0; %detector center offset along the 's' direction in pixels relative to the tube rotation center
 offset_t = -nt/2; % detector center offset along the 't' direction in pixels relative to the tube rotation center
 d_objbottom_det = 0;%in cm,the distance from the bottom of the object to the center detector.
%=======================================
%User defines the recon volume geometry:
%=======================================
%%treat the x-ray tube rotation center as the origin of the 3D coordinate
%%system ( see the coordinate system sketch in the instruction document)
%x: posive direction points toward right, 
%y: positive direction points toward the nipple
%z: positive direction points toward the x-ray source
% voxel size (drx, dry, drz) in cm, 
% dimensions (nrx, nry, nrz) and 
% FOV center offsets (offset_x, offset_y, offset_z): in pixels relative to the rotation center.
% For example, if the coordinates of the FOV center is (xctr, yctr, zctr) relative tothe rotation center,
% then offset_x=-xctr, offset_y=-yctr and offset_z=-zctr.
nrx=504;
nry=156;
drx=0.04;  
dry=drx; 
drz=0.08; 
nrz=64; 
offset_x = 0; %in pixels
offset_y = -nry/2;% in pixels. 0 for full cone, -nry/2 for half cone
zfov = nrz*drz;
offset_z = (dod - (zfov/2 + d_objbottom_det-drz/2))/drz; %in pixels: offset of the volume ctr to the rotation ctr in the z direction; 
                                
 
 %===================
 %Reconstruction
 %===================

 %Generate the system matrix
                          
 igr = image_geom('nx', nrx, 'ny',nry, 'nz', nrz, 'dx',drx, 'dz', drz,...
       'offset_y', offset_y,'offset_z', offset_z,'down', 1); 
 
 
 btg = bt_geom('arc', 'ns', ns, 'nt', nt, 'na', na, ...
		'ds', ds, ...%'dt', dv, ... defautly dt = -ds;
		'down', 1, ...
        'orbit', orbit,...
        'offset_s', 0, ...  
		'offset_t', offset_t, ...
  		'dso', dso, 'dod', dod, 'dfs',inf);  
    
Gtr = Gtomo_syn(btg,igr);

%FBP reconstruction
disp 'FBP'
xfbp = fbp_dbt(Gtr,btg,igr, g,'hann75');

% SART reconstruction
xbp = BP(Gtr, g); %initialization for SART
disp 'SART'
tic
[xartt,costart] = SART_dbt(Gtr,g,xbp,2,0.5);
disp 'SART time '
toc
% ML recosntruction
disp 'ML'
tic
[xmlt,costml] = ML_dbt(Gtr,proj,xbp,I0,3,2);
disp 'ML time'
toc

disp 'Recon completed';


figure('Name','dbr_recon_circular_example');
imagesc(xartt(end:-1:1,:,30)), daspect([1 1 1]), colormap(gray)
title 'Slice 30 (lesion focal plane) of SART reconstruction'
colorbar;

save fbp_cir.mat xfbp;
save sart_cir.mat xartt;
save ml_cir.mat xmlt; 

