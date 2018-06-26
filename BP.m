function img = BP(Gt,proj)
%function img = BP(Gt,proj)
%A simple back-projection reconstruction for DBT.
%
%Inputs:
%   Gt: DBT system projector operator, created using "Gtomo_syn.m".
%   proj: DBT projection views in the form of line integral, a 3D array of 
%       size ns x nt x na, where ns and nt are the projection image dimensions
%       and na is the number of views.
%
%Output:
%   img:  the back-projection reconstructed volume, a 3D array.
%
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018

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


nview =length(Gt);
arg=Gt{1}.arg;

nx=arg.ig.nx;
ny=arg.ig.ny;
nz=arg.ig.nz;
dz=arg.ig.dz;
nmask=sum(arg.ig.mask(:));

ns=arg.cg.ns;%arg.nn(1);
nt=arg.cg.nt;%arg.nn(2);
nd=arg.nd;

img=zeros(nmask,1);
%denom = zeros(1,nmask);
y=zeros(ns,nt);
for i=1:nview
    p=proj(:,:,i);
    G=Gt{i};
%     arg=G.arg;
%     f3d_mex('init', arg.sys_type, uint8(arg.mask), ...
%         int32(arg.nx), int32(arg.ny), int32(arg.nz), ...
%       	 int32(arg.nthread), int32(arg.chat));
     l=sum(G');
     l=reshape(l,ns,nt);
     lmask=logical(l>abs(dz)/10); %Mask out too small l value to avoid blowup in y:
                     %some l can be very small. When it is divided from p, 
                     %the y value can be overamplified at the boundary(the 
                     %boundary value could be 10 or 20 times larger than the
                     %normal object value.
      y=zeros(size(p));
      y(lmask)=p(lmask)./l(lmask);
      clear lmask l
%     y(isnan(y))=0;
     imgi = G'*y(:);
     clear y
     denomi=sum(G);
    
     imgj=imgi./denomi(:);
     clear denomi imgi
     imgj(isnan(imgj))=0;
     img=img+imgj;
     clear imgj
end
%img=embed(img./denom(:),arg.mask);
img=embed(img/nview,arg.ig.mask);
img=permute(img,[1 3 2]);