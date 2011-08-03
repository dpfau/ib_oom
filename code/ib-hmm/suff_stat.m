function [gamma gam_tot xi] = suff_stat(seq,alpha,beta,c,O,T)

g = alpha.*beta;
gamma = zeros(max(seq),size(g,1));
for i = 1:max(seq)
    gamma(i,:) = sum(g(:,seq==i),2)';
end

gam_tot = sum(g(:,2:end),2);

A = join_oom(O,T);
xi = zeros(size(alpha,1),size(alpha,1),max(seq));
for i = 1:max(seq)
    foo = tprod(alpha(:,seq(1:end-1)==i),[1 2],1./c(seq(1:end-1)==i),2);
    bar = beta(:,2:end);
    xi(:,:,i) = tprod(foo,[2 -1],bar(:,seq(1:end-1)==i),[1 -1]).*A(:,:,i);
end