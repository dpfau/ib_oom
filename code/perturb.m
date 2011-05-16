function oom2 = perturb(oom,dx)

oom2 = oom + max(0,dx*randn(size(oom)));
oom2 = oom2./repmat(sum(sum(oom2,3)),[size(oom,1) 1 size(oom,3)]);