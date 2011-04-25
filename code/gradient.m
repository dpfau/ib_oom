function grad = gradient(seq,oom,z_0,b)
% grad(i,j,x) = gradient of log(p(seq)) wrt oom(i,j,x)
% If z_t is the vector of states of an OOM at time t
% grad_z(k,i,j,x) is the derivative of z_t(k) w.r.t. oom(i,j,x),
% where oom(:,:,x) is the observable operator for symbol x, scaled by
% p(seq) so that it does not vanish as t grows large.

nz = size(z_0,1);
nx = size(oom,3);

if numel(oom) ~= nz*nz*nx
    disp('Array of operators is not properly sized.')
end

z = z_0;
q = z_0;
grad_z = zeros(nz,nz,nz,nx);
grad_q = zeros(nz,nz,nz,nx);
log_p = 0;
for i = 1:length(seq)
    grad_1 = tprod(oom, [1 -1 4], grad_z, [-1 2 3 4]);
    grad_1(:,:,:,seq(i)) = grad_1(:,:,:,seq(i)) + tprod(eye(nz), [1 2], z, 3);
    z_1 = oom(:,:,seq(i))*z;
    p = sum(z_1);
    
    grad_z = grad_1/p;
    grad_q = grad_z + grad_q;
    z = z_1/p;
    q = q + z;
    log_p = log_p + log(p);
end

grad = squeeze(sum(grad_q));
%grad = tprod(grad_q, [-1 1 2 3], 1+log(q), -1) + b*squeeze(sum(grad_z));