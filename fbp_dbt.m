function img=fbp_dbt(Gt,btg,ig,proj,window)
% function img=fbp_dbt(Gt,btg,ig,proj,window);
% Filtered backp rojection reconstruction for circular DBT
%
% Inputs:
%   Gt: DBT system projector operator, created using "Gtomo_syn()".
%   btg: A structure variable contains the tomosynthesis system parameters,
%       created using "bt_geom.m'.
%   ig: image geometry structure contains the parameters that describe the 
%       reconstructed volume created using "image_geom.m"in Fessler's IRT
%       package.  
%   proj: DBT projection views in the form of line integral, a 3D array of 
%       size ns x nt x na, where ns and nt are the projection image dimensions
%       and na is the number of views.
%   window: reconstruction filter type, can be 'ramp','hann','hann50',
%       'hann75','hann80' etc, see "fbp2_window.m" in IRT for possible options. 
%
%Output:
%   img:  the FBP reconstructed volume, a 3D array.
%
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018
% Reference: Zeng EtAl, "Evaluating the sensitivity of the optimization of 
% acquisition geometry to the choice of reconstruction algorithm in digital
% breast tomosynthesis through a simulation study", PMB, vol. 60(3), 2015

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
if(nargin<5)
    window = 'hann';
end

orbit = btg.orbit+abs((btg.ad(2)-btg.ad(1)));
nview =btg.na;
if(mod(nview,2))
    ctrview = (nview+1)/2;
else
    error('Not implemented for odd # of views!')
end
arg = Gt{ctrview}.arg;
over=1;
ds = btg.ds/over;
dt = btg.dt/over;
halforbit = btg.orbit/2*pi/180;
orbit_start = btg.ad(1);
dsd = btg.dsd;
dod = btg.dod;
dso = btg.dso;

cbg_det_wid = dod*tan(halforbit) + btg.s(end)/cos(halforbit);
cbg_det_depth = btg.t(end)*(dsd+btg.s(end)*tan(halforbit))/dsd;
ns = ceil(cbg_det_wid/ds)*2*over;
nt = ceil(cbg_det_depth/dt)*2*over;

cbg = ct_geom('fan', 'ns', ns, 'nt', nt, 'na', nview, ...
		'ds', ds, 'dt', dt, ...
        'orbit',orbit, ...
        'orbit_start', orbit_start, ...
		'offset_s', 0, ... % quarter detector
		'offset_t', 0.0, ...
		'dsd', dsd, 'dod', dod, 'dfs', inf);

ang = btg.ad * pi/180 ;
for i=1:nview
    asin=sin(ang(i));
    acos=cos(ang(i));
    xs = -dso*asin;
    ys = 0;
    zs = dso*acos;
    [u,v] = meshgrid(cbg.s,cbg.t);
    
    xd = u*acos + dod*asin;
    yd = v;
    zd = u*asin - dod*acos;

    k = (-dod-zd)./(zs-zd);
    
    xi = xd + k.*(xs-xd);
    yi = yd + k.*(ys-yd);
    
   % [newxi, newyi]=meshgrid(xi,yi);
    [oldxi,oldyi]=meshgrid(arg.cg.s, arg.cg.t);
    p=interp2(oldxi,oldyi, (proj(:,:,i))', xi,yi,'linear');
    pcb(:,:,i)=p';
end
pcb(isnan(pcb))=0;
img = feldkamp(cbg, Gt{ctrview}.arg.ig, pcb, 'use_mex', 0,'window',window);
img = permute(img,[1 3 2]);
img = img *180/btg.orbit; %scale  

%for test
% ellt=[0 ig.offset_z*dz 0 8 5 6 0 0 0.1];
% projt = ellipsoid_proj(cbg, ellt);
% imgt = feldkamp(cbg, Gt{ctrview}.arg.ig, projt, 'use_mex', 0);    
%     
%     
    