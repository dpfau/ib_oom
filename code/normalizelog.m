function y = normalizelog(x)

expx = exp(x - ones(size(x,1),1)*max(x));
y = expx./(ones(size(x,1),1)*sum(expx));