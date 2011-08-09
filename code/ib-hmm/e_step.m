function [o,t] = e_step(g,gt,xi,O,T,Z,beta)

eps = 1e-6;
err = Inf;

n = size(T,1);
k = size(T,3);

o = O;
t = T;
z = Z;

while err > eps
    %% Find best Lagrange multipliers
    A = [zeros(n*k,1), kron(ones(k,1),eye(n)), zeros(n*k), reshape(t,n,n*k)'.*repmat(z,k,n); ...
        zeros(n^2*k,1+n), kron(eye(n*k),ones(n,1)), kron(o(:).*repmat(z,k,1),eye(n)); ...
        ones(n,1), zeros(n,n*(k+1)), tprod(t,[2 1 -1],o,[1 -1]) - eye(n)];
    
    b = [beta * g(:) ./ o(:); ...
        xi(:) ./ t(:); ...
        - gt(:) ./ z(:)];
    b(isnan(b)) = 0;
    
    l = pinv(A)*b;
    err_l = norm(A*l - b);
    
    l_O = l(2:n+1);
    l_T = l(n+1+(1:n*k));
    l_z = l(2+n*(k+1):end);
    
    %% Find best stationary distribution
    bar = [];
    for i = 1:k
        bar = [bar; diag(l_z'*t(:,:,i))];
    end
    foo = [];
    for i = 1:k
        foo = [foo; diag(o(:,i))];
    end
    A = [tprod(o(:),1,bar,[1 2]); ...
        tprod(t(:),1,kron(foo,l_z),[1 2]); ...
        diag(-l(1) + l_z - tprod(t,[1 2 -1],o,[2 -1])*l_z); ...
        tprod(t,[1 2 -1],o,[2 -1]) - eye(n)];
    
    b = [beta * g(:) - kron(ones(k,1),l_O) .* o(:); ...
        xi(:) - kron(l_T,ones(n,1)) .* t(:); ...
        gt(:); ...
        zeros(n,1)];
    
    z = pinv(A)*b;
    err_z = norm(A*z - b);
    
    %% Find best emission distribution
    foo = [];
    for i = 1:k
        foo = [foo; l_O + z.*(t(:,:,i)'*l_z)];
    end
    bar = [];
    for i = 1:k
        bar = [bar, diag(l_z'*t(:,:,i))];
    end
    A = [diag(foo); ...
        tprod(t(:),1,kron(diag(repmat(z,k,1)),l_z),[1 2]); ...
        tprod(z,1,bar,[1 2]);...
        repmat(eye(n),1,k)];
    b = [beta * g(:); ...
        xi(:) - kron(l_T,ones(n,1)) .* t(:); ...
        -gt(:) + (l_z - l(1)) .* z; ...
        ones(n,1)];
    
    o = pinv(A)*b;
    err_o = norm(A*o - b);
    o = reshape(o,n,k);
    
    %% Find best transition distribution
    foo = [];
    for i = 1:k
        foo = [foo, kron(diag(o(:,i)),l_z')];
    end
    A = [kron(diag(repmat(z,k,1)),l_z'); ...
        diag(kron(l_T,ones(n,1)) + repmat(l_z,n*k,1) .* kron(repmat(z,k,1) .* o(:),ones(n,1))); ...
        foo;...
        kron(eye(n*k),ones(1,n))];
    b = [beta * g(:) - kron(ones(k,1),l_O) .* o(:); ...
        xi(:); ...
        -gt(:) + (l_z - l(1)) .* z; ...
        ones(n*k,1)];
    
    t = pinv(A)*b;
    err_t = norm(A*t - b);
    t = reshape(t,n,n,k);
    
    %% Display the error
    fprintf('Error - L: %d, Z: %d, O: %d, T: %d\n',err_l,err_z,err_o,err_t);
end