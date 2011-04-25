function [x fval] = easymin(seq,n,b)

nx = max(seq);
z_0 = zeros(n,1);
z_0(1) = 1;

f = @(x) objective(seq,reshape(x,[n n nx]),z_0,b);
Aeq = kron(repmat(eye(n),1,nx),ones(1,n));

x0 = rand(n,n,nx);
x0 = x0./repmat(sum(sum(x0,3)),[n 1 nx]);

[x fval] = fmincon(f,x0(:),[],[],Aeq,ones(n,1),zeros(n*n*nx,1),ones(n*n*nx,1));