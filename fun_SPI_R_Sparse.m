function [im_r, totaliter] = fun_SPI_R_Sparse(patterns, measurements, para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, Oct 21, 2017
% Contact: lihengbian@gmail.com
% This function implements the singel pixel imaging reconstruction using
% the total variation compressive sensing method (lagrange multiplier method).
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

mu_0 = 1e-3;
rho = 1.05;
plot_flag=0;

tol=1e-2; % accuracy
min_iter = 30;
x0 = ones(row * col,1);
if exist('para','var')
    if isfield(para,'tol')
        tol = para.tol; % accuracy
    end
    if isfield(para,'min_iter')
        min_iter = para.min_iter;
    end
    if isfield(para,'x0')
        x0 = para.x0;
    end
    if isfield(para,'mu')
        mu_0 = para.mu;
    end
end

%%
% CS reconstruction
H = generate_H_DCT(row,col);
[x_tv, totaliter]=m_tv(H,P,measurements,x0,mu_0,rho,tol,min_iter,plot_flag);

im_r = reshape(x_tv,[row,col]);
end



%% subfuntion generate the transform matrix H (DCT basis)
function [H] = generate_H_DCT(M,N)
% generate DCT (discrete cosine transform) operator matrix

mm = 0:M-1;
mm = repmat(mm,[1,N]); % m coordinate
nn = zeros(1,M*N);
for temp = 0:N-1
    nn(1,temp*M+1:temp*M+M) = temp;
end

H = zeros(M*N,M*N);
k = 0;
for p = 0:M-1
    for q = 0:N-1
        if p == 0
            ap = 1/sqrt(M);
        else
            ap = sqrt(2/M);
        end
        if q == 0
            aq = 1/sqrt(N);
        else
            aq = sqrt(2/N);
        end
        k = k + 1;
        H(k,:) = ap .* aq .* cos(pi*(2*mm+1)*p/2/M).*cos(pi*(2*nn+1)*q/2/N);
    end
end

end



%% subfunction m_tv
function [x, totaliter] = m_tv(H,P,y,x0,mu_0,rho,tol,min_iter,plot_flag)
% solve the following model using the Lagrange multiplier method
% min ||Hx||_1
% s.t.Px=y;

x=x0;
w=H*x;
r=zeros(size(w));
s=zeros(size(y));
object_kp=norm(P*x-y);
o=object_kp;
mu=mu_0;
mu_bar=mu_0;
converged=false;
iter=0;
max_iter=300;

% % pinv_mat=pinv(H'*H+P'*P);

% speed up (by Liheng, 2016-06-01)
temp = H'*H+P'*P;
pinv_mat=temp\eye(size(temp));
pinv_matH = pinv_mat*H';
pinv_matP = pinv_mat*P';

if plot_flag==1
    figure;
    norm_w_hx=norm(w-H*x);
    norm_px_y=norm(P*x-y);
    norm_x=[norm(x)];
    norm_w=[norm(w)];
    norm_r=[norm(r)];
    norm_s=[norm(s)];
end

while (~converged) && (iter<=max_iter)
    iter=iter+1;
    
    % update w
    tmp=H*x-1/mu*r;
    flag1=tmp>1/mu;
    flag2=tmp<-1/mu;
    w=(tmp-1/mu).*flag1+(tmp+1/mu).*flag2;
    
    % update x
    x = pinv_matH*(w+1/mu*r) + pinv_matP*(y-1/mu*s);
    
    % update r & s
    r=r+mu*(w-H*x); 
    s=s+mu*(P*x-y);
     
    % update mu
    mu=min(mu_bar,mu*rho);   
    
    % stop criterion
    object_k=object_kp;
    object_kp=norm(P*x-y);
    stop_c1=object_k-object_kp;    
    if (stop_c1<tol && stop_c1>=0) && iter>min_iter
        converged=true;
% %         disp('iteration is converged');
    elseif iter>max_iter
        converged=true;
% %         disp('iter number reached maximun number');
    end
    
    if mod(iter,10) == 0
        fprintf(['Sparse iteration ' num2str(iter) ', the error is ' num2str(stop_c1) '\n']);
    end
    
% %     if plot_flag==1
% % %             norm_x=[norm_x,norm(x)];
% % %             norm_w=[norm_w,norm(w)];
% % %             norm_r=[norm_r,norm(r)];
% % %             norm_s=[norm_s,norm(s)];
% % %             norm_w_hx=[norm_w_hx norm(w-H*x)]
% % %             norm_px_y=[norm_px_y norm(P*x-y)]
% %             o=[o,norm(H*x,1)];
% % %             subplot(241);plot(norm_w_hx);title('w-hx');
% % %             subplot(242);plot(norm_px_y);title('px-y');
% % %             subplot(243);plot(norm_x);title('x');
% % %             subplot(244);plot(norm_w);title(['w' sprintf(' mu %1.3f',mu)]);
% % %             subplot(245);plot(norm_r);title('norm r');
% % %             subplot(246);plot(norm_s);title('norm s');
% %             subplot(121);plot(o);title('object value');
% %             subplot(122);imshow(reshape(x,[N1,N2]));title('rst img');
% %             pause(0.01);
% %     end
end

totaliter = iter;
fprintf(['Sparse total iterations ' num2str(totaliter) ', the error is ' num2str(stop_c1) '\n']);

end