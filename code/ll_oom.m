 function ll = ll_oom(data,oom,p)

c = ones(length(data),1); % observation probabilities
if nargin > 2
    z = p;
else
    z = stat(oom);
end

z = diag(sum(oom(:,:,data(1))))*z;
c(1)   = sum(z);
z = z/c(1);
for t = 2:length(data)
    z    = oom(:,:,data(t))*z;
    c(t) = sum(z);
    z    = z/c(t);
end
ll = sum(log(c));
