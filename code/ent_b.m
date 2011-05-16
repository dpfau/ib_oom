function ent = ent_b(oom)
% for a given oom (or oom representation of an hmm), approximate H[B] by
% generating data and using the nearest neighbors entropy estimate

zs = zeros(size(oom,1),1e5);

z = stat(oom);
p_xz = squeeze(sum(oom))';
for i = 1:1e5
    prob = p_xz*z;
    cusum = cumsum(prob);
    x = find(cusum >= cusum(end)*rand,1);
    z = oom(:,:,x)*z;
    z = z/sum(z);
    zs(:,i) = z;
end

ent = nn_ent(zs');