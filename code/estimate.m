function [p q] = estimate(seq,oom,z)
% Estimates p(x_t+1|x_t,z_t-1) and p(x_t|z_t-1)from data

nz = size(oom,1);
nx = size(oom,3);

if nargin == 2
    z = stat(oom);
end
p = zeros(nx,nx,nz);
q = zeros(nx,nz);

for i = 1:length(seq)-1
    p(seq(i+1),seq(i),:) = p(seq(i+1),seq(i),:) + permute(z, [2 3 1]);
    q(seq(i),:) = q(seq(i),:) + z';
    z = oom(:,:,seq(i))*z;
    z = z/sum(z);
end
q(seq(end),:) = q(seq(end),:) + z';

p = tprod(p, [1 2 3], 1./squeeze(sum(p)), [2 3],'n'); % turn the joint into the conditional
q = q./(ones(nx,1)*sum(q));