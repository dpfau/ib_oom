function z_0 = stat(oom)
% returns the stationary distribution of the OOM

T = sum(oom,3); % transition matrix for latent states
if nnz(T<0) == 0
    dz = 1;
    z_0 = rand(size(oom,1));
    z_0 = z_0/sum(z_0);
    while dz > 1e-9
        z_1 = T*z_0;
        dz = max(abs(z_1-z_0));
        z_0 = z_1;
    end
else
    [v,d] = eig(T);
    z_0 = v(:,find(diag(d)>1-1e-5 & diag(d)<1+1e-5,1));
    z_0 = z_0/sum(z_0);
end