function [im_r, totaliter] = fun_SPI_R_CGD(patterns, measurements, para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 22, 2017
% Contact: lihengbian@gmail.com
% This function implements the singel-pixel imaging reconstruction using the conjugate gradient descent method .
% If this code offers any help, please cite the publication:
% Liheng Bian, Jinli Suo, Qionghai Dai, and Feng Chen. 'Experimental comparison of single-pixel imaging algorithms'.

% Inputs:
% patterns: illumination patterns (pixels * pixels * pattern numbers)
% measurements: single pixel measurements (vector)

% Outputs:
% im_r: reconstructed image (pixels * pixels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row, col, m] = size(patterns);
P = reshape(patterns, [row*col, m]);
P = P'; % each row represents a pattern

tol=1e-2; % accuracy
min_iter = 30;
x = ones(row * col,1);
if exist('para','var')
    if isfield(para,'tol')
        tol = para.tol; % accuracy
    end
    if isfield(para,'min_iter')
        min_iter = para.min_iter;
    end
% %     if isfield(para,'x0')
% %         x = para.x0;
% %     end
end

%%
A = P'*P; % row*col * row*col
b = P'*measurements; % row*col

% solve Ax = b
N = row*col;

rr = b - A*x; % error
r = rr;
d = rr; % gradient

for k=0:N-1
    a = (norm(rr)^2)/(d'*A*d);   % step
    x = x + a * d;   % update
    
    rr = b-A*x;    % error
    rrr = abs(norm(r)-norm(rr));
    
    B=(norm(rr)^2)/(norm(r)^2);
    d=rr+B*d;
    
    r=rr;
    
    if mod(k+1,100) == 0
        fprintf(['CGD iteration ' num2str(k+1) ', the error is ' num2str(rrr) '\n']);
    end
    
    if ((rrr<=tol)||(k==N-1)) && k>min_iter
        break;
    end
end

im_r = reshape(x,[row,col]);

totaliter = k;
fprintf(['CGD total iterations ' num2str(totaliter) ', the error is ' num2str(rrr) '\n']);

end
