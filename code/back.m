function beta = back(seq,oom)
% Just the backward pass of the forward-backward algorithm

beta = ones(size(oom,1),length(seq));
for t = length(seq)-1:-1:1
    beta(:,t) = beta(:,t+1)'*oom(:,:,seq(t+1));
    beta(:,t) = beta(:,t)/sum(beta(:,t));
end