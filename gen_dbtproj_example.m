%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example code to simulate DBT projection views.
% 
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%attenuation properties of breast tissue: Jonn and
%Yaffe-1987-pmb-v32-p678-Table1
mu_fg=0.802;%0.378 @30keV;  0.802 @20keV;
mu_adp=0.456;%0.264 @30keV; 0.456 @20keV
mu_carci=0.844;% carcinoma 0.392 @30keV; 0.844 @20keV;
mu_ca=1.2; %calcification


load example_data/breast_phantom0.mat; %the variable name is 'x'
 %downsample to run faster
 down=2;
 x = downsample3(x,down);
%===============================    
%Define the object geometry
%===============================
[nx,ny,nz]=size(x);%phantom dimensions
dx=0.02*down; dy=0.02*down; dz=0.02*down; %in cm, phantom pixel sizes
xfov = dx*nx;%20;
yfov = dy*ny;
zfov = dz*nz;

offset_y = -ny/2;% offset of the object ctr in pixels for the y-dimension; 0 for full cone, -ny/2 for half cone
d_objbottom_det = 0; %in cm,the distance from the bottom of the object to the center detector.
                     %Value "0" means the object is places right on the
                     %detector.

%============================
%Define the scanner geometry
%============================
% for arc trajectory
 dso = 60.5; %in cm: dist. from the source to the rotation center
 dod = 4.5;  %in cm: dist. from the rotation center to the detector
 dsd = dso+dod; %in cm: dist. from source to the detector
 
 
 orbit =30;  %angular span
 na = 7; %number of projection views
 ds = dx; % in cm; detector pixel size
 dt= dx;
 %calculate the length and width of the detector so it is large enough to cover the
 %projection views from the most oblique angles.
 %costheta=cos(orbit/2*pi/180); sintheta=sin(orbit/2*pi/180);
 %sfov = ((dso*costheta+dod)*(xfov/2+dso*sintheta)/(dso*costheta+offset_z*dz-zfov/2) - dso*sintheta)*2;
 %tfov = yfov*(dso*costheta+dod)/(dso*costheta+dod-zfov);
 ns = 615;%ceil(sfov/ds);
 nt = 170;%ceil(tfov/dt);

 offset_s = 0; %detector center offset along the 's' direction in pixels relative to the tube rotation center
 offset_t = -nt/2; %detector center offset along the 't' direction in pixels relative to the tube rotation center
 offset_z = (dod - (zfov/2 + d_objbottom_det - dz/2 ))/dz; %in pixels, offset of the object ctr to the rotation ctr in the z direction;

 %==============================
 %Create DBT projection views
 %==============================
 btg = bt_geom('arc', 'ns', ns, 'nt', nt, 'na', na, ...
		'ds', ds, ...%'dt', dv, ... defautly dt = -ds;
		'down', 1, ...
        'orbit', orbit,...
		'offset_s', 0, ... % quarter detector 
		'offset_t', offset_t, ...
   		'dso', dso, 'dod', dod, 'dfs',inf);  

ig = image_geom('nx', nx, 'ny',ny, 'nz', nz, 'dx',dx, 'dz', dz,...
       'offset_y', offset_y,'offset_z', offset_z,  'down', 1); 

Gt = Gtomo_syn(btg,ig); %generate system Fatrix

%add a spherical lesion to the phantom
if(1)
    rx=dx*8; ry=dy*8; rz=dz*4;%lesion radius
    %define the geometric properties of a sphere
    %ell=[xctr, yctr, zctr, rx, ry, rz, alpha, beta, attenuation];
    ell=[ig.x(252) ig.y(78) ig.z(60) rx ry rz 0 0 1];      
    lesion = ellipsoid_im(ig,ell); %generate the lesion volumes
    x(lesion==1) = 3*mu_ca; %assign the attnuation value to lesion voxels.
                            %This lesion is bright so will be highly
                            %visuable in the recosntructed DBT volume.
end

nview=length(Gt);
g=zeros(btg.ns,btg.nt,nview);
for i=1:nview
   tic   
   g(:,:,i)=Gt{i}*permute(x,[1 3 2]);   
   toc
end

%===========================
%add Poisson noise
%===========================
g_noi=g;

I0=3*10^5/na; %distriburte the entire dose evenly to each proejction view. 
proj=I0*exp(-g);
proj_noi=proj;
%add poisson noise to the projections
if(1)
    proj_noi=poissrnd(proj);%poisson(proj);%
    proj_noi(find(proj_noi==0))=1; 
    
    g_noi=log(I0)-log(proj_noi); % convert back to line integrals 
    g_noi(g_noi<0)=0;    
end

if(1)
    save proj_noi.mat proj_noi I0;
    save g_noi.mat g_noi;
end


