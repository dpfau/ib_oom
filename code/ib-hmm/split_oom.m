function [O,T] = split_oom(oom)

O = squeeze(sum(oom,1));
T = tprod(oom,[1 2 3],1./O,[2 3]);