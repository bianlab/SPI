function [im_r, totaliter] = fun_SPI_R_Poisson(patterns, measurements, para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 22, 2017
% Contact: lihengbian@gmail.com
% This function implements the singel pixel imaging reconstruction using the Poisson maximum likelihood method.
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
alpha = 0.1;
beta = 0.5;

A = P;
b = measurements;

k = 0;
r = norm(b - A*x);

% % figure; hold on;

while 1
    k = k+1;
    g = grad(A,x,b);
    
   % initialize step size
    temp = x./g; % keep the step not so large, so that x>0
    temp(temp<=0) = 1;
    step = min(1,min(temp));
    
    % search for the step size
    ktemp = 1;
    while funL(A,x-step*g,b) > funL(A,x,b) - alpha*step*g'*g
        ktemp = ktemp + 1;
        step = step * beta;
    end
    x = x - step * g;
    
    rlast = r;
    r = norm(b - A*x);
    
    % show results
    if mod(k,100) == 0
        fprintf(['Poisson iteration ' num2str(k) ', the error is ' num2str(rlast-r) '\n']);
        im_r = reshape(x,[row,col]);
% %         imshow(im_r,[],'InitialMagnification',1000); title(['Poisson iteration ' num2str(k)]);
    end
    
    % if reaches convergence
    if (rlast-r<tol || k>3*row*col) && k>min_iter
        break        
    end
end

im_r = reshape(x,[row,col]);

totaliter = k;
fprintf(['Poisson total iterations ' num2str(totaliter) ', the error is ' num2str(rlast-r) '\n']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = funL(A,x,b)
    temp = ones(1,size(A,1));
    ax = A*x;
    L = temp*(ax-b.*abs(log(ax)));
end

function g = grad(A,x,b)
% %     g = A'*((A*x-b)./(A*x));
    g = A'*(1 - b./(A*x));
end
