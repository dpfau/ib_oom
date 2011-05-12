function ll = score(dat,oom,z_0)

if nargin == 2
    z = stat(oom);
else
    z = z_0;
end

ll = cell(size(dat));
for t = 1:length(dat)
    seq = dat{t};
    ll{t} = zeros(size(seq));
    for i = 1:length(seq)
        z = oom(:,:,seq(i))*z;
        p = sum(z);
        z = z/p;
        ll{t}(i) = p;
    end
end