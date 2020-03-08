function [im_r, totaliter] = fun_SPI_R_GD(patterns, measurements, para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 22, 2017
% Contact: lihengbian@gmail.com
% This function implements the singel pixel imaging reconstruction using the gradient descent method .
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
    if isfield(para,'x0')
        x = para.x0;
    end
end

%%
A = P;
b = measurements;
AA = A'*A;

r = b - A*x;
normr = norm(r);

k = 0;
while 1
    k = k+1;
    g = grad(A,x,b,AA);
    step = -(g'*A'*r)/(g'*(AA)*g);
    x = x - step * g;
    r = b - A*x; 
    
    normrlast = normr;
    normr = norm(r);
    
    if mod(k,100) == 0
        fprintf(['GD iteration ' num2str(k) ', the error is ' num2str(normrlast - normr) '\n']);
    end
    if (normrlast - normr<tol || k>3*row*col) && k>min_iter
        break        
    end
end

im_r = reshape(x,[row,col]);

totaliter = k;
fprintf(['GD total iterations ' num2str(totaliter) ', the error is ' num2str(normrlast - normr) '\n']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = grad(A,x,b,AA)
% %     g = 2*A'*(A*x-b);
    g = 2*(AA*x-A'*b);
end

