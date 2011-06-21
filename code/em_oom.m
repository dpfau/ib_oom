function oom = em_oom(seq,test,n,oom0)

if nargin < 4
    oom = rand(n,n,max(cellfun(@max,seq)));
    oom = oom./repmat(sum(sum(oom,3)),[size(oom,1) 1 size(oom,3)]);
else
    oom = oom0;
end
    
p = rand(n,1);
p = p/sum(p);

l = -1;
ll = -1;
while (ll == -1) || ((l - ll)/sum(cellfun(@length,seq)) > 1e-4)
    ll = l;
    [oom,p,l,i,j] = pass(seq,oom,p);
    fprintf('Train LL: %d, Test LL: %d, I[X;Z]: %d, I[Z;Y]: %d, H[Z]: %d\n',...
        l/sum(cellfun(@length,seq)),...
        sum(cellfun(@(x) ll_oom(x,oom),test))/sum(cellfun(@length,test)),...
        i,...
        j,...
        entropy(stat(oom)));
end

function [oom,p,ll,ixz,izy] = pass(seq,OOM,P)
% One pass of the EM algorithm for an OOM

c = cell(size(seq));
oom = zeros(size(OOM));
p = zeros(size(P));
ixz = 0;
izy = 0;
for i = 1:length(seq)
    [a,b,g,c{i},xi] = for_back_oom(seq{i},OOM,P);
    oom = oom + xi;
    p = p + g(:,1);
    ixz = ixz - sum(entropy(a));
    b2 = b .* (stat(OOM)*ones(1,size(b,2)));
    b2 = (ones(size(b,1),1)*(1./sum(b2))) .* b2;
    izy = izy - sum(entropy(b2));
end
oom = oom./repmat(sum(sum(oom,3)),[size(oom,1) 1 size(oom,3)]);
p = p/sum(p);
ll = sum( cellfun( @(x) sum( log( x ) ), c ) );

ixz = ixz/sum(cellfun(@length,seq));
izy = izy/sum(cellfun(@length,seq));

ixz = ixz + entropy(stat(OOM));
izy = izy + entropy(stat(OOM));