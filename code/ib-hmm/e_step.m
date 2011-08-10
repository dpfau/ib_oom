function [o,t] = e_step(g,gt,xi,beta,O,T,Z)

len = numel(O) + numel(T) + numel(Z);
n = size(T,1);
k = size(T,3);

assert( size(T,2) == n );
assert( size(O,1) == n );
assert( size(O,2) == k );
assert( size(Z,1) == n );
assert( size(Z,2) == 1 );

x = fmincon( @(x) [beta * g(:); xi(:); -gt(:)]'*log(x), ...
    [O(:); T(:); Z(:)], ...
    [], ...
    [], ...
    [kron(ones(1,k),eye(n)), zeros(n,n*(n*k+1)); ...
     zeros(n*k), kron(eye(n*k),ones(1,n)), zeros(n*k,n); ...
     zeros(1,n*k*(n+1)), ones(1,n)], ...
    ones(n+n*k+1,1), ...
    zeros(len,1), ...
    ones(len,1), ...
    @(x) stat_con(x,n,k), ...
    optimset('Algorithm','interior-point'));

o = reshape(x(1:n*k),n,k);
t = reshape(x(n*k+(1:n^2*k)),n,n,k);

function [C Ceq] = stat_con(x,n,k)
% constraint that Z is stationary distribution of 

o = reshape(x(1:n*k),n,k);
t = reshape(x(n*k+(1:n^2*k)),n,n,k);
z = x((n+1)*n*k+(1:n));

C = 0;
Ceq = (tprod(t,[1 2 -1],o,[2 -1]) - eye(n)) * z;