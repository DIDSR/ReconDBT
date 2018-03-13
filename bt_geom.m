 function st = bt_geom(type, varargin)
% function st = bt_geom(type, varargin)
% Create the DBT geometry structure that descrives the system geometry.
% It is modified based on Jeff Fessler's ct_geom() in the IRT package.
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018
%----------------------------------------------------
%ct_geom():
% Create the "ct geometry" structure that describes the sampling
% characteristics of a tomosyntheis system.
% --for 'lin' (lineartomo): coordinate center is the fulcrum
%
% in:
%%	type	 'lin' (linear tomo)
%
% options for all geometries
%%	'orbit_start'	starting angle for breastCT ('ct')	
%%                  center position of the source orbit for linear tomo
%%                  ('lin')
%%	'orbit'			distance that the source moves
%					
%	'down'			down-sampling factor, for testing
%	'units'			string to print distance units (default: 'mm')
%
% options for fan-beam
%%	'ns'			# of horizontal samples
%%	'nt'			# of vertical samples
%%	'na' | 'nbeta'		# of samples
%%	'ds'			horizontal sample spacing (default: 1)
%%	'dt'			vertical sample spacing (default: -ds)
%	'offset_s'		unitless (default: 0)
%			(relative to centerline between two central channels).
%			use 0.25 or 1.25 for "quarter-detector offset"
%	'offset_t'		unitless (default: 0)
%
% options for helical
%	'zshifts' [na]		user-specified source z positions.
%				usually this is empty (default) in which case
%				it is computed internally from "ztrans/z_dir"
%				or "pitch" (preferred).  has units!
%	'pitch'			bed_travel_per_rotation / axial_fov. default: 0.
%				(unitless. can be negative. usually near 1.0)
%	'ztrans'		modified total translation of the source in z.
%				ztrans = dz_source * na. default:0 (axial).
%				Modified because dz_source = pitch / na_per_rot.
%					Not recommended!
%	'z_dir'			default:1 i.e. source moves from negative to
%					positive z quadrant. -1 for opposite.
%					Not recommended!
%	'offset_z'		z-location of the 1+floor(na/2) detector,
%				in number of dz_source units (unitless).
%				Default 0.  Often +/- 0.5, or table reference.
%
%	fan beam distances:
%%	'dsd' | 'dis_src_det'	Dist from the source plane to the detector
%%	plane
%%	'dso' | 'dis_src_iso'	Dist from the fulcrum to the source plane
%%	'dod' | 'dis_iso_det'	Dist from the fulcrum to the detector plane
%%	'dfs' | 'dis_foc_src'	default: inf (flat detector)
%                           use 0 for 3rd generation CT arc
%
% out:
%	st	(struct)	initialized structure
%
% methods:
%	st.shape(sino)		reshape sinograms that are columns into 3d array
%	st.s			s sample locations
%	st.t			t sample locations
%	st.ws			(ns-1)/2 + st.offset_s
%	st.wt			(nt-1)/2 + st.offset_t
%%	st.ad			source angles in degrees (for breast CT) or 
%%                  displacement in mm (for linear tomo)
%	st.ar			source angles in radians
%%  st.dd           detector displacement in mm
%	st.dim			dimensions: [st.ns st.nt st.na]
%	st.downsample(down)	reduce sampling by integer factor
%	st.ones			ones(ns,nt,na)
%	st.zeros		zeros(ns,nt,na)
%	st.rmax			max radius within FOV
%	st.shape(sino(:))	reshape to [ns,nt,na,?]
%	st.unitv(is,it,ia)	unit 'vector' with single nonzero element
%	st.plot([ig])		show geometry
%
%	trick: you can make orbit=0 and orbit_start = column vector (length na)
%	if you need nonuniformly spaced projection view angles.
%
% Copyright 2006-1-18, Jeff Fessler, The University of Michigan
%-----------------------------------------------------------------

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                       Legal Disclaimer
% This software and documentation (the "Software") were developed at the 
% Food and Drug Administration (FDA) by employees of the Federal Government 
% in the course of their official duties. Pursuant to Title 17, Section 105
% of the United States Code, this work is not subject to copyright protection
% and is in the public domain. Permission is hereby granted, free of charge, 
% to any person obtaining a copy of the Software, to deal in the Software 
% without restriction, including without limitation the rights to use, copy, 
% modify, merge, publish, distribute, sublicense, or sell copies of the 
% Software or derivatives, and to permit persons to whom the Software is 
% furnished to do so. FDA assumes no responsibility whatsoever for use by 
% other parties of the Software, its source code, documentation or compiled 
% executables, and makes no guarantees, expressed or implied, about its 
% quality, reliability, or any other characteristic. Further, use of this 
% code in no way implies endorsement by the FDA or confers any advantage in 
% regulatory decisions. Although this software can be redistributed and/or 
% modified freely, we ask that any derivative works bear some notice that 
% they are derived from it, and any modified versions bear some notice that 
% they have been modified.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin == 1 && streq(type, 'test'), ct_geom_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

if streq(type, 'ge1') % special case
	st = ct_geom_ge1(type, varargin{:});
return
end

% defaults
st.type = type;
st.ns = [];
st.nt = [];
st.na = [];
st.down = 1;
st.orbit_start = 0;
st.pitch = 0;
st.ztrans = 0;
st.offset_z = 0;
st.z_dir = 1;
st.units = 'mm';
st.user_zshifts = [];

if (streq(type, 'arc') | streq(type, 'lin'))
	st = ct_geom_fan(st, varargin{:});
%elseif streq(type, 'par')
%	st = ct_geom_par(st, varargin{:});
%elseif streq(type, 'moj')
%	st = ct_geom_moj(st, varargin{:});
else
	error('unknown sinotype %s', type)
end

if isempty(st.na), st.na = 2 * floor(st.ns * pi/2 / 2); end

ct_geom_zshift_check(st)

meth = { ...
	's', @ct_geom_s, '()'; ...
	't', @ct_geom_t, '()'; ...
	'ws', @ct_geom_ws, '()'; ...
	'wt', @ct_geom_wt, '()'; ...
	'ad', @ct_geom_ad, '()'; ...
	'dd', @ct_geom_ar, '()'; ... 
    'zshifts', @ct_geom_zshifts, '()'; ...
	'downsample', @ct_geom_downsample, '()'; ...
	'dim', @ct_geom_dim, '()'; ...
	'ones', @ct_geom_ones, '()'; ...angular sample
	'rmax', @ct_geom_rmax, '()'; ...
	'unitv', @ct_geom_unitv, '() | (is,it,ia)'; ...
	'zeros', @ct_geom_zeros, '()'; ...
	'shape', @ct_geom_shape, '()'; ...
	'plot', @ct_geom_plot, '() | (ig)';
	};

st = strum(st, meth);

if st.down ~= 1
	down = st.down; st.down = 1; % trick
	st = st.downsample(down);
end


% ct_geom_dim()
function dim = ct_geom_dim(st)
dim = [st.ns st.nt st.na];
if isempty(st.ns) || isempty(st.nt) || isempty(st.na)
	error 'dim requested without ns,nt,na'
end


% ct_geom_ones()
% sinogram of all ones
function out = ct_geom_ones(st)
out = ones(st.dim);


% ct_geom_unitv()
% sinogram with a single ray
function out = ct_geom_unitv(st, is, it, ia)
out = st.zeros;
if ~isvar('is') || isempty(is)
	is = floor(st.ns/2 + 1);
	it = floor(st.nt/2 + 1);
	ia = 1;
end
out(is,it,ia) = 1;


% ct_geom_zeros()
% sinogram of all zeros
function out = ct_geom_zeros(st)
out = zeros(st.dim);


% ct_geom_rmax()
% max radius within fov
function rmax = ct_geom_rmax(st)
smax = max(abs(st.s));
if streq(st.type, 'fan')
	if isinf(st.dso) % parallel
		rmax = smax;angular sample
	elseif st.dfs == 0 % arc
		rmax = st.dso * sin(smax / st.dsd);
	elseif isinf(st.dfs) % flat
		rmax = st.dso * sin(atan(smax / st.dsd));
	else
		error 'unknown case'
	end
end


% ct_geom_ws()
% 'middle' sample position
function ws = ct_geom_ws(st)
ws = (st.ns-1)/2 + st.offset_s;


% ct_geom_wt()
% 'middle' sample position
function wt = ct_geom_wt(st)
wt = (st.nt-1)/2 + st.offset_t;


% ct_geom_s()
% sample locations ('radial')

function s = ct_geom_s(st, varargin)
s = st.ds * ([0:st.ns-1]' - st.ws);
if length(varargin)
	s = s(varargin{:});
end

% ct_geom_t()
% sample locations ('radial')
function t = ct_geom_t(st, varargin)
t = st.dt * ([0:st.nt-1]' - st.wt);
if length(varargin)
	t = t(varargin{:});
end


% ct_geom_ad()
% source sampling locations (degrees for 'ct' and mm for 'lin')
function ang = ct_geom_ad(st, varargin)
if streq(st.type,'arc')
    %ang = [0:st.na-1]'/st.na * st.orbit + st.orbit_start;
    ang = [-(st.na-1)/2: (st.na-1)/2] * st.orbit/(st.na-1) + st.orbit_start;
    ang = ang(varargin{:});
elseif streq(st.type,'lin')
    ang = [-(st.na-1)/2: (st.na-1)/2] * st.orbit/(st.na-1) + st.orbit_start;
    ang = ang(varargin{:});
else 
    disp 'unknow ct type!'
    return;
end

% ct_geom_dd()
%function ang = ct_geom_dd(st, varargin)
% if streq(st.type,'arc')
%     ang = st.dod*tan(ct_geom_ad(st)*pi/180); % deg2rad(ct_geom_ad(st));
%     ang = ang(varargin{:});
% elseif streq(st.type,'lin')   
%     ang = ct_geom_ad(st) *(-st.dod/st.dso);
%     ang = ang(varargin{:});
% else
%    disp 'unknow ct type!'
%     return;
%end

% ct_geom_ar()
% angular sample locations (radians)
function ang = ct_geom_ar(st, varargin)
ang = deg2rad(ct_geom_ad(st));
ang = ang(varargin{:});



% ct_geom_shape()
% reshape into sinogram array
function sino = ct_geom_shape(st, sino)
sino = reshape(sino, st.ns, st.nt, st.na, []);


% ct_geom_downsample()
% down-sample (for testing)
function st = ct_geom_downsample(st, down)
st.down = st.down * down;

if ~isempty(st.user_zshifts)
	st.user_zshifts = st.user_zshifts(1:down:st.na);
end

st.ns = 2 * round(st.ns / down / 2); % keep it even
st.nt = 2 * round(st.nt / down / 2); % keep it even
st.na = length([1:down:st.na]);

if streq(st.type, 'fan')
	st.ds = st.ds * down;
	st.dt = st.dt * down;
else
	error(['unknown sinotype ' type])
end


% ct_geom_zshift_check()
function ct_geom_zshift_check(st)
if ~isempty(st.user_zshifts) + (st.pitch ~= 0) + (st.ztrans ~= 0) > 1
 fail('only one of "pitch" (recommended) or ztrans/z_dir or user zshifts can be given')
end
if ~isempty(st.user_zshifts) && length(st.user_zshifts) ~= st.na
	fail('user zshifts size mismatch')
end

%
% ct_geom_zshifts()
%
% source z locations
function zshifts = ct_geom_zshifts(st, varargin)

if ~isempty(st.user_zshifts) % user-provided
	zshifts = st.user_zshifts;
	zshifts = zshifts(varargin{:});
return
end

if st.pitch % by jeff, tested in ellipsoid_proj.m
	na_per_360 = st.na * 360 / st.orbit;
	D = st.dso / (st.dso + st.dod) * st.nt * st.dt; % collimaty = ((0:ig.ny-1) - ig.ny/2) * ig.dy;
	dz_source = st.pitch * D / na_per_360; % source shift per view
	zshifts = (st.offset_z + [0:st.na-1]' - floor(st.na/2)) * dz_source;

else % by somesh
	dz_source = st.ztrans / st.na; % possibly 0 if axial
	if st.z_dir == 1
		zshifts = (st.offset_z + [0:st.na-1]' - floor(st.na/2)) * dz_source;
	else
		zshifts = (st.offset_z + [st.na-1:-1:0]' - (st.na - floor(st.na/2) - 1)) * dz_source; % fix: not tested
	end
end
zshifts = zshifts(varargin{:});


%
% ct_geom_fan()
%
function st = ct_geom_fan(st, varargin);

% defaults
st.orbit = 360; % [degrees]
st.ds		= 1;
st.dt		= [];
st.offset_s	= 0;
st.offset_t	= 0;

st.dsd = [];	% dis_src_det
st.dso = [];	% dis_src_iso
st.dod = [];	% dis_iso_det
st.dfs = 0;	% dis_foc_src (3rd gen CT)

subs = { ...
	'src_det_dis', 'dsd';
	'dis_src_det', 'dsd';
	'dis_src_iso', 'dso';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs';
	'nbeta', 'na';
	'zshifts', 'user_zshifts';
	};
st = vararg_pair(st, varargin, 'subs', subs);
if isempty(st.dt), st.dt = -st.ds; end

% work out distances
if (~isempty(st.dsd) && isinf(st.dsd)) ...
|| (~isempty(st.dso) && isinf(st.dso)) % handle parallel-beam case gracefully
	st.dsd = inf; st.dso = inf; st.dod = 1;
end
if isempty(st.dsd) + isempty(st.dso) + isempty(st.dod) > 1
	error 'must provide at least two of dsd, dso, dod'
end
if isempty(st.dsd), st.dsd = st.dso + st.dod; end
if isempty(st.dso), st.dso = st.dsd - st.dod; end
if isempty(st.dod), st.dod = st.dsd - st.dso; end
if st.dso + st.dod ~= st.dsd
	error 'bad fan-beam distances'
end


%
% ct_geom_plot2()
% picture of 2D source position / detector geometry
%
function ct_geom_plot2(st, ig)
if ~streq(st.type, 'fan'), error 'only fan done', end
x0 = 0;
y0 = st.dso;
t = linspace(0,2*pi,1001);
switch st.dfs
case 0
	gam = st.s / st.dsd; % 3rd gen: equiangular
case inf
	gam = atan(st.s / st.dsd); % flat
otherwise
	error 'not done'
end
xds = st.dsd * sin(gam);
yds = st.dso - st.dsd * cos(gam);
rot = deg2rad(st.orbit_start);
rot = [cos(rot) sin(rot); -sin(rot) cos(rot)];
p0 = rot * [x0; y0];
pd = rot * [xds'; yds'];
rfov = st.dso * sin(max(abs(gam)));

plot(	0, 0, '.', ...
	p0(1), p0(2), 's', ...
	[pd(1,1) p0(1) pd(1,end)], [pd(2,1) p0(2) pd(2,end)], '-', ...
	st.dso * cos(t), st.dso * sin(t), '--', ... % source circle
	rfov * cos(t), rfov * sin(t), ':', ...
	pd(1,:), pd(2,:), 'o')
if isvar('ig') && ~isempty(ig)
	hold on
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	plot([xmax xmin xmin xmax xmax], [ymax ymax ymin ymin ymax], 'g-')
	hold off
end
title(sprintf('fov = %g', rfov))
axis square, zoom on


%
% ct_geom_plot3()
% picture of 3D helical CT geometry
%
function ct_geom_plot3(st, ig)
if ~streq(st.type, 'fan'), error 'only fan done', end

t1 = -st.dso * sin(st.ar); % src x pos
t2 = st.dso * cos(st.ar); % src r pos
t3 = st.zshifts; % src z pos
plot3(-st.dso * sin(st.ar), st.dso * cos(st.ar), st.zshifts, 'y.') % trajectory
view(-110, 22)
hold on
text(t1(1), t2(1), t3(1), '1')
center = 1 + floor(st.na/2);
text(t1(center), t2(center), t3(center), sprintf('%d',center))
text(t1(end), t2(end), t3(end), sprintf('%d',st.na))
ar
% axes
plot3([-100*ceil(max(abs(t1))/100) 100*ceil(max(abs(t1))/100)], [0 0], [0 0])
text(100*ceil(max(abs(t1))/100),0, 'x')
plot3([0 0], [-100*ceil(max(abs(t2))/100) 100*ceil(max(abs(t2))/100)], [0 0])
text(0, 100*ceil(max(abs(t2))/100), 'y')
plot3([0 0], [0 0], [-10*ceil(max(abs(t3))/10) 10*ceil(max(abs(t3))/10)])
text(0, 0, 10*ceil(max(abs(t3))/10), 'z')

grid on
xlabel(sprintf('x (%s)', st.units))
ylabel(sprintf('y (%s)', st.units))
zlabel(sprintf('z (%s)', st.units))

if isvar('ig') && ~isempty(ig) && isfield(ig.meth,'z')
	hold on
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	zmin = min(ig.z); zmax = max(ig.z);
	plot3(	[xmin xmax xmax xmin xmin ], ...
		[ymin ymin ymax ymax ymin ], ...
		[zmin zmin zmin zmin zmin ], ...
		'g-')
	plot3(	[xmin xmax xmax xmin xmin], ...
		[ymin ymin ymax ymax ymin], ...ar
		[zmax zmax zmax zmax zmax], ...
		'g-')
	plot3(xmin*[1 1], ymin*[1 1], [zmin zmax], 'g-')
	plot3(xmax*[1 1], ymin*[1 1], [zmin zmax], 'g-')
	plot3(xmin*[1 1], ymax*[1 1], [zmin zmax], 'g-')
	plot3(xmax*[1 1], ymax*[1 1], [zmin zmax], 'g-')
end

% detector array for source at first position
cen_loc = [t1(1) t2(1) t3(1)];
unit_vec = [-cen_loc(1) -cen_loc(2) 0] / sqrt(cen_loc(1)^2+cen_loc(2)^2);
% diametrically opposite point of the source on the detector
det_cen_loc = cen_loc + st.dsd * unit_vec;
plot3(det_cen_loc(1),det_cen_loc(2),det_cen_loc(3),'rs')
theta = st.s/st.dsd;
det_x = cen_loc(1) + st.dsd * (unit_vec(1)*cos(theta)-unit_vec(2)*sin(theta));
det_y = cen_loc(2) + st.dsd * (unit_vec(1)*sin(theta)+unit_vec(2)* ...
				cos(theta));
st_t = st.t;
for tt=1:length(st_t)
	plot3(det_x, det_y, cen_loc(3)+st_t(tt)*ones(size(theta)), 'r.')
end

% detector array for source at middle position
cen_loc = [t1(center) t2(center) t3(center)];
unit_vec = [-cen_loc(1) -cen_loc(2) 0] / sqrt(cen_loc(1)^2+cen_loc(2)^2);

% diametrically opposite point of the source on the detector
det_cen_loc = cen_loc + st.dsd * unit_vec;
plot3(det_cen_loc(1),det_cen_loc(2),det_cen_loc(3),'rs')
theta = st.s/st.dsd;
det_x = cen_loc(1) + st.dsd * (unit_vec(1)*cos(theta)-unit_vec(2)*sin(theta));
det_y = cen_loc(2) + st.dsd * (unit_vec(1)*sin(theta)+unit_vec(2)* ...
				cos(theta));
st_t = st.t;
for tt=1:length(st_t)
	plot3(det_x, det_y, cen_loc(3)+st_t(tt)*ones(size(theta)), 'c.')
end

% detector array for source at final position
cen_loc = [t1(end) t2(end) t3(end)];
unit_vec = [-cen_loc(1) -cen_loc(2) 0] / sqrt(cen_loc(1)^2+cen_loc(2)^2);
% diametrically opposite point of the source on the detector
det_cen_loc = cen_loc + st.dsd * unit_vec;
plot3(det_cen_loc(1),det_cen_loc(2),det_cen_loc(3),'rs')
theta = st.s/st.dsd;
det_x = cen_loc(1) + st.dsd * (unit_vec(1)*cos(theta)-unit_vec(2)*sin(theta));
det_y = cen_loc(2) + st.dsd * (unit_vec(1)*sin(theta)+unit_vec(2)* ...
				cos(theta));
st_t = st.t;
for tt=1:length(st_t)
	plot3(det_x, det_y, cen_loc(3)+st_t(tt)*ones(size(theta)), 'm.')
end


%
% ct_geom_plot()
%
function out = ct_geom_plot(st, varargin)
im clf
if all(st.zshifts == 0) % 2D or axial case
	ct_geom_plot2(st, varargin{:});
else
	ct_geom_plot3(st, varargin{:});
end
if nargout, out = []; end


%
% ct_geom_test()
%
function ct_geom_test
% axial cone-beam test
cg = ct_geom('fan', 'ns', 888, 'nt', 64, 'na', 984, ...
	'offset_s', 1.25, ...
	'dsd', 949, 'dod', 408);
cg.ad(2);
cg.downsample(2);
cg.rmax;
cg.ws;
cg.s(cg.ns/2+1);
cg.plot;
prompt

% test user zshifts
cg = ct_geom('fan', 'dsd', 949, 'dod', 408, ...
	'ns', 888, 'ds', 1.0239, 'offset_s', 1.25, ...
	'nt', 8, 'dt', 1.0964, ...
	'na', 2*984, 'orbit', 2*360, ...
	'zshifts', []);
cg.zshifts(1:2:7);

% helical cone-beam test
cg = ct_geom('fan', 'dsd', 949, 'dod', 408, ...
	'ns', 888, 'ds', 1.0239, 'offset_s', 1.25, ...
	'nt', 8, 'dt', 1.0964, ...
	'na', 2*984, 'orbit', 2*360, ...
	'pitch', 1.0);
%	'ztrans', 16);

cg.zshifts(1:2:7);
cg.downsample(2);
cg.s(cg.ns/2+1);
cg.plot;
