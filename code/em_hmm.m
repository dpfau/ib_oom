function [t,o] = em_hmm(seq,test,n,t0,o0)

if nargin < 4
    t = rand(n);
    t = t./(ones(size(t,1),1)*sum(t));
else
    t = t0;
end
   
if nargin < 5
    o = rand( max( cellfun( @max, seq ) ), n );
    o = o./(ones(size(o,1),1)*sum(o));
else
    o = o0;
end
    
p = rand(n,1);
p = p/sum(p);

l = -1;
ll = -1;
while (ll == -1) || ((l - ll) > 1e-12)
    ll = l;
    [t,o,p,l] = pass(seq,t,o,p);
    fprintf('Training LL: %d, Testing LL: %d\n',l/sum(cellfun(@length,seq)),sum(cellfun(@(x) ll_hmm(x,t,o,p),test))/sum(cellfun(@length,test)));
end

function [t,o,p,ll] = pass(seq,T,O,P)
% One pass of the EM algorithm for a hidden markov model

g = cell(size(seq));
c = cell(size(seq));

t = zeros(size(T));
o = zeros(size(O));
p = zeros(size(P));
for i = 1:length(seq)
    [~,~,g{i},c{i},xi] = for_back(seq{i},T,O,P);
    t = t + xi;
    for j = 1:size(O,1)
        o(j,:) = o(j,:) + sum(g{i}(:,seq{i}==j),2)';
    end
    p = p + g{i}(:,1);
end
t = t./(ones(size(t,1),1)*sum(t));
o = o./(ones(size(o,1),1)*sum(o));
p = p/sum(p);
ll = sum( cellfun( @(x) sum( log( x ) ), c ) );
