function oom = rand_oom(nx,nz)

oom = betarnd(.1, 10, nz, nz, nx);
norm = sum(sum(oom,3));
oom = oom./repmat(norm,[nz 1 nx]);