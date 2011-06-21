function ib_hmm(seq,n,b,t0,o0)

if nargin < 4
    t = normalize(rand(n));
else
    t = t0;
end
   
if nargin < 5
    o = normalize(rand(max(seq),n));
else
    o = o0;
end

while 1
    oom = zeros(n,n,size(o,1));
    for i = 1:size(o,1)
        oom(:,:,i) = diag(o(i,:))*t;
    end
    [t,o] = pass(seq,oom,b);
    subplot(2,1,1); imagesc(t)
    subplot(2,1,2); imagesc(o); drawnow
end

function [t,o] = pass(seq,oom,b)

beta = back(seq,oom);
alpha = normalize(diag(stat(oom))*beta.^-b);
gamma = normalize(alpha.*beta);

t = zeros(size(oom,1),size(oom,2));
for i = 1:length(seq)-1
    t = t + ( (beta(:,i+1)*alpha(:,i)') .* oom(:,:,seq(i+1)) )/sum(oom(:,:,seq(i))*alpha(:,i));
end

o = zeros(size(oom,3),size(oom,1));
for i = 1:max(seq)
    o(i,:) = o(i,:) + sum(gamma(:,seq==i),2)';
end

t = normalize(t);
o = normalize(o);