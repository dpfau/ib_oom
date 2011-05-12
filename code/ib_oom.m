function p_z1_x1__z0 = ib_oom(seq,seq_,nz,b)
% Learning Hidden Markov Models by Information Bottleneck
% David Pfau, 2011

% The convention used is that double underscores represent conditioning, so
% p_z1__x1_z0 = p(z_t+1 | x_t+1, z_t) whereas
% p_z1_x1__z0 = p(z_t+1, x_t+1 | z_t), which are the observable operators.

if strcmp(class(seq),'double')
    nx = max(seq);
else
    nx = 1;
    for t = 1:length(seq)
        nx = max(nx,max(seq{t}));
    end
end

p_z1_x1__z0 = rand_oom(nx,nz);
F = 100;
while 1
    p_z0 = stat(p_z1_x1__z0);
    [p_x2__x1_z0 p_x1__z0] = estimate(seq, p_z1_x1__z0, p_z0);
    p_z1__x1_z0 = tprod(p_z0, 1, exp( -b*d_kl(p_x2__x1_z0, p_x1__z0) ), [1 2 3], 'n');
    
    Z = squeeze(sum(p_z1__x1_z0));
    F1 = -tprod(log(Z), [-1 -2], p_x1__z0, [-2 -1], 'n');
    if abs(F-F1) < 1e-9, break; end
    F = F1;
    
    p_z1__x1_z0 = p_z1__x1_z0./tprod(ones(nz,1), 1, Z, [2 3], 'n'); % normalize
    p_z1_x1__z0 = tprod(p_z1__x1_z0, [1 2 3], p_x1__z0, [3 2]);
    l = log_p(seq,p_z1_x1__z0);
    fprintf('Free Energy: %4.2d, Data Log Likelihood: %4.2d, Test Log Loss: %4.2d\n',F,l,log_p(seq_,p_z1_x1__z0));
end

function d = d_kl(p,q)
% d(i,j,k) = D_kl[p(:,k,j)||q(:,i)]
d = zeros(size(q,2),size(p,3),size(p,2));
for i = 1:size(q,2)
    for j = 1:size(p,3)
        for k = 1:size(p,2)
            log_pq = log(p(:,k,j)./q(:,i));
            log_pq(isinf(log_pq) | isnan(log_pq)) = 0;
            d(i,j,k) = p(:,k,j)'*log_pq;
        end
    end
end