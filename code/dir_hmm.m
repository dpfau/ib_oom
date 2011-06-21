function [hmm O T] = dir_hmm(n,m,a,b)
% Returns a random HMM with 'n' observables and 'm' latent states, in the form of an observable operator model.
% Each observation probability p( x_t | z_t ) is sampled from a uniform Dirichlet with 'a' total pseudocounts
% Each transition probability p( z_{t+t} | z_t ) is sampled from a uniform Dirichlet with 'b' total pseudocounts

T = zeros(m);
O = zeros(n,m);

for i = 1:m
    O(:,i) = gamrnd(a/n*ones(n,1),1);
    O(:,i) = O(:,i)/sum(O(:,i));
    
    T(:,i) = gamrnd(b/m*ones(m,1),1);
    T(:,i) = T(:,i)/sum(T(:,i));
end

addpath '/Users/davidpfau/Documents/MATLAB/tprod'
hmm = tprod(T,[1 2],O,[3 2],'n');