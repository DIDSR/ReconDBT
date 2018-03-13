%Run setpath_DBT to include the mfiles in Fessler's IRT reconstruction packages that are used in ReconDBT
 
%You should change the path of "irt" and "reconDBT" to the actual path
%where you stored the packages.
irt='/raida/rpz/matlab/irt/'; %The path for Fessler's IRT package 
reconDBT = '/raida/rpz/matlab/dbt/ReconDBT/'; % The path for the reconDBT package 

addpath(reconDBT);
addpath([irt 'mex/v7']);

list = {...
'fbp', ...		% FBP (filtered backprojection) code
'systems', ...		% system "matrices"
'utilities', ...	% various utility functions
};

for ii=1:numel(list)
	tmp = [irt list{ii}];
	if exist(tmp, 'dir'), addpath(tmp), end
end

if strcmp([irt 'fbp' filesep 'cbct_back.m'], which('cbct_back'))
	disp('Path setup for irt appears to have succeeded.')
	clear list ii irt tmp
else
	disp('Path setup for irt may have failed.')
end


disp 'Completed!'
