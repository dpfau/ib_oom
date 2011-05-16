function h = nn_ent( x )
% Nearest neighbors estimate of entropy for a continuous distribution
% x - NxK matrix of N K-dimensional vectors, i.e. each row is a point

addpath '/Users/davidpfau/Documents/Code/kdtree1.2/kdtree'
tree = kdtree_build( x );

h = 0;
for i = 1:size(x,1)
    nn = kdtree_k_nearest_neighbors( tree, x(i,:)', 2 );
    h = h + 1/2 * max( -1000, log( ( x(i,:) - x(nn(1),:) ) * ( x(i,:) - x(nn(1),:) )' ) ); % max is to prevent numerical underflow
end
h = h/size(x,1) + log(size(x,1)) + log(2) - psi(1);
kdtree_delete( tree );