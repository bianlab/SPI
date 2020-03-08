function [im_r] = fun_SPI_R_DGI(patterns, measurements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 25, 2016
% Contact: lihengbian@gmail.com
% This function implements the singel pixel imaging reconstruction using the differential ghost imaging (DGI) method proposed in 
% [1] W. Gong and S. Han. A method to improve the visibility of ghost images obtained by thermal light. Phys. Lett. A, 374(8):1005-1008, 2010.
% [2] F. Ferri, D. Magatti, L.A. Lugiato, and A. Gatti. Differential ghost imaging. Phys. Rev. Lett., 104(25):253603, 2010.

% Inputs:
% patterns: illumination patterns (pixels * pixels * pattern numbers)
% measurements: single pixel measurements (vector)

% Outputs:
% im_r: reconstructed image (pixels * pixels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
B_aver  = 0;
SI_aver = 0;
R_aver = 0;
RI_aver = 0;
iter = 0;

%% DGI reconstruction
for j = 1:size(patterns,3)
    
    pattern = patterns(:,:,j);
    
    iter = iter + 1;

    B_r     = measurements(j);

    % Differential ghost imaging (DGI)
    SI_aver = (SI_aver * (iter -1) + pattern * B_r)/iter;
    B_aver  = (B_aver * (iter -1) + B_r)/iter;
    R_aver = (R_aver * (iter -1) + sum(sum(pattern)))/iter;
    RI_aver = (RI_aver * (iter -1) + sum(sum(pattern))*pattern)/iter;

    im_r = SI_aver - B_aver / R_aver * RI_aver;
end
end