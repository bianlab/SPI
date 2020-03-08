function [im_r, totaliter] = fun_SPI_R_AP(patterns, measurements, para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 25, 2016
% Contact: lihengbian@gmail.com
% This function implements the singel pixel imaging reconstruction using the alternating projection method proposed in 
% "Kaikai Guo, Shaowei Jiang, and Guoan Zheng, Multilayer fluorescence imaging on a single-pixel detector, Biomedical Optics Express, 7, 7, 2425 (2016).".

% Inputs:
% patterns: illumination patterns (pixels * pixels * pattern numbers)
% measurements: single pixel measurements (vector)

% Outputs:
% im_r: reconstructed image  (pixels * pixels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row, col, m] = size(patterns);
P = reshape(patterns, [row*col, m]);
P = P'; % each row represents a pattern
A = P;
b = measurements;

tol=1e-2; % default accuracy
min_iter = 30; % default minimum iteration
x = ones(row * col,1);
if exist('para','var')
    if isfield(para,'tol')
        tol = para.tol; % accuracy
    end
    if isfield(para,'min_iter')
        min_iter = para.min_iter;
    end
    if isfield(para,'x0')
        x = para.x0;
    end
end

%%
im_r = reshape(x,[row, col]);
r = b - A*im_r(:); 

for j = 1 : 3*row*col
    % begin iteration
    for i = 1 : size(patterns,3)
        temp_im = im_r .* patterns(:,:,i);
        temp_im1 = temp_im;    
        temp_im = temp_im - mean(mean(temp_im)) + measurements(i)/numel(temp_im);
        
        im_r = abs(im_r + (patterns(:,:,i))./((max(max(patterns(:,:,i)))).^2) .* (temp_im - temp_im1));
    end;
    
    rlast = r;
    r = b - A*im_r(:);
    
    if mod(j,100) == 0
        fprintf(['AP iteration ' num2str(j) ', the error is ' num2str((norm(rlast) - norm(r))) '\n']);
    end
    if (abs(norm(rlast) - norm(r))<tol || j>3*row*col) && j>min_iter
        break        
    end
    
end

totaliter = j;
fprintf(['AP total iterations ' num2str(totaliter) ', the error is ' num2str((norm(rlast) - norm(r))) '\n']);

end