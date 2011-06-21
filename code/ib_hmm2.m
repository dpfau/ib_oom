function ib_hmm2(seq,n,b,a)

if nargin < 4
    oom = rand(n,n,max(seq));
    oom = oom./repmat(sum(sum(oom,3)),[size(oom,1) 1 size(oom,3)]);
else
    oom = a;
end

while 1
    [oom l] = pass(seq,oom,b);
    for i = 1:size(oom,3)
        subplot(size(oom,3),1,i); 
        imagesc(oom(:,:,i)); 
        drawnow
    end
    fprintf('LL: %d, Lagrangian: %d\n',ll_oom(seq,oom)/length(seq),l);
end

function [oom_,l] = pass(seq,oom,b)

beta = back(seq,oom);
alpha = normalizelog(log(stat(oom))*ones(1,size(beta,2)) - b * log(beta));

oom_ = zeros(size(oom));
for i = 1:length(seq)-1
    oom_(:,:,seq(i+1)) = oom_(:,:,seq(i+1)) + ( (beta(:,i+1)*alpha(:,i)') .* oom(:,:,seq(i+1)) )/sum(oom(:,:,seq(i))*alpha(:,i));
end
oom_ = oom_./repmat(sum(sum(oom_,3)),[size(oom_,1) 1 size(oom_,3)]);
c = beta .* (stat(oom_)*ones(1,size(beta,2)));
c = (ones(size(c,1),1)*(1./sum(c))) .* c;
l = (1-b)*entropy(stat(oom_)) + mean(entropy(alpha)) - b*mean(entropy(c));