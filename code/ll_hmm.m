 function ll = ll_hmm(data,T,O,p)

if nargin > 3
    for j = 1:size(O,1)
        A{j} = sparse(diag(O(j,:)))*T;
    end
else
    for j = 1:size(T,3)
        A{j} = T(:,:,j);
    end
    p = O;
end
c = ones(length(data),1); % observation probabilities
if nargin > 2
    z = p;
else
    z = stat(T);
end

z = diag(sum(A{data(1)}))*z;
c(1)   = sum(z);
z = z/c(1);
for t = 2:length(data)
    z = A{data(t)}*z;
    c(t)   = sum(z);
    z = z/c(t);
end
ll = sum(log(c));