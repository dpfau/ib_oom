function [curve p_zx p_z p_yz] = ib(p_x,p_yx,k,b0,b1,dT)
% Implements the iterative update for the information bottleneck
% Solves for p(z|x), p(z), p(y|z) given p(x) and p(y|x) 
% with k codewords and inverse temperature b
% From Tishby, Pereira and Bialek 1999
% David Pfau, 2011

sx = size(p_yx,2);
Hx  = -p_x'*log(p_x);
Ixy = I(p_x,p_yx);

p_zx = rand(k,sx); % initialize signal-to-codeword mapping randomly
obj = Inf;
b = b0;
curve = [];
while 1
    p_zx = p_zx./(ones(k,1)*sum(p_zx)); % normalize
    p_z  = p_zx*p_x;
    p_yz = p_yx*(p_zx.*((1./p_z)*p_x'))';
    imagesc(p_zx); drawnow
    
    Iyz = I(p_z,p_yz);
    Izx = I(p_x,p_zx);
    if abs(obj - (Izx - b*Iyz)) < 1e-8
        if b <= b1
            break; 
        else
            curve = [curve; b, Iyz/Ixy, Izx/Hx];
            fprintf('Beta: %f -> %f\n',b,1/((1/b) + dT));
            b = 1/((1/b) + dT);
            
            dp_zx = rand(k,sx);
            dp_zx = dp_zx./(ones(k,1)*sum(dp_zx));
            l = .5;
            
            p_zx = l*p_zx + (1-l)*dp_zx; % convex combination of current solution and random probability vectors;
            p_z  = p_zx*p_x;
            p_yz = p_yx*(p_zx.*((1./p_z)*p_x'))';
            
            Iyz = I(p_z,p_yz);
            Izx = I(p_x,p_zx);
        end
    end
    obj  = Izx - b*Iyz;
    
%     subplot(3,1,1); imagesc(p_zx); title('p(z|x)')
%     subplot(3,1,2); imagesc(p_zx.*((1./p_z)*p_x'));  title('p(x|z)')
%     subplot(3,1,3); imagesc(p_yz); title('p(y|z)'); drawnow
        
    p_zx = p_z*ones(1,sx).*exp(-b*d_kl(p_yx,p_yz));
    F = -log(sum(p_zx))*p_x; % the free energy
    fprintf('Free energy: %2.4f, Objective: %2.4f, I[Z;Y]/I[X;Y]: %1.4f, I[Z;X]/H[X]: %1.4f\n',F,obj,Iyz/Ixy,Izx/Hx);
end

function d = d_kl(p,q)
% d(i,j) = D_kl[p(:,j)||q(:,i)]
d = zeros(size(q,2),size(p,2));
for i = 1:size(q,2)
    for j = 1:size(p,2)
        d(i,j) = p(:,j)'*log(p(:,j)./q(:,i));
    end
end

function i = I(p_x,p_yx)
% i = I[X;Y] = H[Y] - H[Y|X]
p_y = p_yx*p_x;
log_p = log(p_yx./(p_y*ones(1,size(p_yx,2))));
log_p(isnan(log_p) | isinf(log_p)) = 0;
i = sum(p_yx.*log_p)*p_x;