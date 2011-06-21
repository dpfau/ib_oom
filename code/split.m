function [O2 T2] = split(O,T)

O2 = [O O];
rnd = rand(size(T,1),2*size(T,2));
T2 = cat(1,rnd .* [T T], (1-rnd) .* [T T]);