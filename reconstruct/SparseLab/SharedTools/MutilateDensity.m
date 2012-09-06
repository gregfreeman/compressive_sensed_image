function Mf = MutilateDensity(grid,f,a,b)
% MultilateDensity -- set density to zero outside (a,b)
%  Mf = MutilateDensity(grid,f,a,b)
Mf = f .* ((a <= grid) & (grid <= b));

%
% Copyright (c) 2006. David Donoho
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
