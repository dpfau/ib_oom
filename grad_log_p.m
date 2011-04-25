function grad = grad_log_p(seq,oom,z_0)
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

dz = zeros(nz,nz,nz,nx);
z  = z_0;
grad = zeros(nz,nz,nx);

for i = 1:length(seq)
    dz1 = tprod(oom, [1 -1 4], dz, [-1 2 3 4]);
    dz1(:,:,:,seq(i)) = dz1(:,:,:,seq(i)) + tprod(eye(nz), [1 2], z, 3);
    z1 = oom(:,:,seq(i))*z;
    p = sum(z1);
    dz = dz1/p - tprod(z1, 1, squeeze(sum(dz1)), [2 3 4])/p^2;
    z = z1/p;
    grad = grad + squeeze(sum(dz1))/p;
end