function [x z] = gen_hmm(O,T,len,p)
% Generates data of length len from HMM with observation matrix O,
% transition matrix T, initial probabilit p.

if nargin < 4
    zt = sample(stat(T));
else
    zt = sample(p);
end

x = zeros(len,1);
z = zeros(len,1);
for i = 1:len
    x(i) = sample(O(:,zt));
    z(i) = zt;
    zt   = sample(T(:,zt));
end