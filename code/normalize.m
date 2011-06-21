function y = normalize(x)

y = x./(ones(size(x,1),1)*sum(x));