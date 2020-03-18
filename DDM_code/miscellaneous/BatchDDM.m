% runs sequence of functions for analysis

% clearvars % -global

fs = ometiffreader;  % read in image stack

tic;Iqtau = DDM_core(fs);toc  % do Fourier Transforms, with timing

fit_Iqtau_colloids          % Do fitting of I(tor(q)) \propto tor(q) over q 

% refit_Iqtau_colloids      % Do plots (but not the fits!) again