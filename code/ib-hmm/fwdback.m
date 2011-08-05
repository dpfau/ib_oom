function [alpha beta c] = fwdback(seq,O,T)

if nargin == 3
    A = join_oom(O,T);
else
    A = O;
end

alpha = zeros(size(A,1),length(seq));
beta  = zeros(size(A,1),length(seq));
c = zeros(length(seq),1);

z = stat(A);
for i = 1:length(seq)
    alpha(:,i) = z;
    z = A(:,:,seq(i))*z;
    c(i) = sum(z);
    z = z/c(i);
end

beta(:,end) = O(:,seq(end))/c(end);
for i = length(seq)-1:-1:1
    beta(:,i) = A(:,:,seq(i))'*beta(:,i+1)/c(i);
end