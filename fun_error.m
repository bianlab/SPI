function [imError] = fun_error(groundtruth, im_r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, Oct 24, 2017
% Contact: lihengbian@gmail.com
% This function calculates the normalized error (RMSE) between groundtruth and im_r.
% If this code offers any help, please cite the publication:
% Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

% Inputs:
% groundtruth (pixels * pixels)
% im_r (pixels * pixels)

% Outputs:
% imError: relative error between two images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imError = norm(groundtruth - im_r,'fro')/sqrt(size(groundtruth,1)*size(groundtruth,2));
imError = imError/mean(mean(groundtruth));

end