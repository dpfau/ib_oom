function [seq log_p] = gen_oom(oom,len,z_0)

if nargin == 2 % initialize to the stationary distribution of the latent states
   z_0 = stat(oom);
end

seq = zeros(len,1);
z = z_0;
log_p = 0;
p_xz = squeeze(sum(oom))';
for i = 1:len
    prob = p_xz*z;
    cusum = cumsum(prob);
    seq(i) = find(cusum >= cusum(end)*rand,1);
    z = oom(:,:,seq(i))*z;
    log_p = log_p + log(sum(z));
    z = z/sum(z);
end