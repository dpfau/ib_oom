function obj = objective(seq,oom,z_0,b)

l = 0; % total log probability of the sequence
z = z_0; % normalize probability of the latent state
q = z_0; % sum of normalized probabilities of latent state

for i = 1:length(seq)
    z = oom(:,:,seq(i))*z;
    p = sum(z);
    z = z/p;
    q = q + z;
    l = l + log(p);
end

obj = l;
%obj = -q'*log(q) - b*l;