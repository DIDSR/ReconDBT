function [img,cost] = SART_dbt(Gt,proj,x0,niter,stepsize,saveiter)
% function [img, cost] = SART_dbt(Gt,proj,x0,niter,stepsize,saveiter)
% SART reconstruction for DBT.
% Inputs:
%   Gt: DBT system projector operator, created using "Gtomo_syn()".
%   proj: DBT projection views in the form of line integral, a 3D array of 
%       size ns x nt x na, where ns and nt are the projection image dimensions
%       and na is the number of views.
%   x0: an initial estimate of the reconstuction volume.
%   niter: number of iterations for SART.
%   stepsize: the stepsize for each update.
%   saveiter: 0 or 1. if "0",only output the reconstruction of the final
%           iteration; if "1", output the reconstruction at all iterations.
%
%Outputs:
%   img: the SART reconstructed volume. 
%       If "saveiter" is 0, then "img" is a 3D array contains the volume 
%       of the final iternation;
%       If "saveiter" is 1, then "img" is 1 4D array contains the volumes 
%       of all the iterations.
%   cost: a niterx1 vector containing the cost function value at each iteration
%
% Author: Rongping Zeng, FDA/CDRH/OSEL/DIDSR, 
% Contact: rongping.zeng@fda.hhs.gov
% Feb. 2018
% Reference: Zeng EtAl, "Evaluating the sensitivity of the optimization of 
% acquisition geometry to the choice of reconstruction algorithm in digital
% breast tomosynthesis through a simulation study", PMB, vol. 60(3), 2015.

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

if(nargin<6)
    saveiter=0;
end
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

if(saveiter)
    img=zeros([size(x0) niter]);
else
    img=zeros([size(x0) 1]);
end

y=zeros(ns,nt);

x=permute(x0,[1 3 2]);
clear x0
for iter=1:niter
    iter
    for i=1:nview
        
        G=Gt{i};
%         arg=G.arg;
%         f3d_mex('init', arg.sys_type, uint8(arg.mask), ...
%             int32(arg.nx), int32(arg.ny), int32(arg.nz), ...
%              int32(arg.nthread), int32(arg.chat));
         l=sum(G');
         l=reshape(l,ns,nt);
         lmask=logical(l>5*dz);%Mask out too small l value to avoid blowup in y:
                     %some l can be very small. When it is divided from p, 
                     %the y value can be overamplified at the boundary(the 
                     %boundary value could be 10 or 20 times larger than the
                     %normal object value.
         pdif=proj(:,:,i)-G*x;
         y(lmask)=pdif(lmask)./l(lmask);
         clear l lmask
    %     y(isnan(y))=0;
         imgi = G'*y(:);
         clear y
         denomi=sum(G);
         imgj=imgi./denomi(:);
         clear imgi
         imgj(isnan(imgj))=0;
         x=x+stepsize*reshape(imgj,size(x)); 
         clear imgj
    end
    x(x<0)=0;
    %calculate the cost
    if(nargout>1)
        mse=0;
        for i=1:nview
            pdif=proj(:,:,i)-Gt{i}*x;
            mse=mse + sum(pdif(:).^2);
        end
        cost(iter)=mse;
    end     
    if (saveiter)
        img(:,:,:,iter)=permute(x,[1 3 2]);
    end
end
%img=embed(img./denom(:),arg.mask);
%img=embed(x,arg.mask);
if(~saveiter)
    img=permute(x,[1 3 2]);
end
