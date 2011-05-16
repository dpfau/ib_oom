function [n,h] = cover(x,dx)
% Finds the number (not necessarily minimal!) of balls of radius dx needed
% to cover all the points in x (counting each column as a point), as well
% as the expected code length of a vector encoded to precision dx
% n - covering number
% h - expected code length, should be equal to H - log(n) (+ const due
% to difference between sphere and cube volume?) where H is entropy of p(x)

n = 0;
h = 0;
len = length(x);
while ~isempty(x)
    pt = ceil(rand*size(x,2));
    nn = sum((x - x(:,pt)*ones(1,size(x,2))).^2) > dx^2; % nearest neighbors
    nnn = nnz(~nn); % number of nearest neighbors
    x = x(:,nn);
    n = n + 1;
    h = h - nnn/len*log(nnn/len);
end
h = h + (n-1)/(2*len); % Miller-Madow bias correction