function [ll p] = log_p(dat,oom,z_0)

if nargin == 2
    z = stat(oom);
else
    z = z_0;
end

l = 0;
len = 0;
p = zeros(length(dat),1);
for t = 1:length(dat)
    seq = dat{t};
    len = len + length(seq);
    for i = 1:length(seq)
        z = oom(:,:,seq(i))*z;
        p(i) = sum(z);
        l = l + log(p(i));
        z = z/p(i);
    end
end

ll = -l/len/log(2);