 function [f b g c xi] = for_back(data,T,O,p)
% Given data, observation matrix O and transition matrix T, returns the
% forward vectors f, backward vectors b, and the marginal vectors g (gamma)
% as well as the probability of observations c

% If we want to pass an OOM insted of HMM, just leave out the third
% argument, and pass the 3D array of linear observable operators as T.
    
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
f = zeros(size(T,1),length(data)); % forward probabilities
b = zeros(size(T,1),length(data)); % backward probabilities
c = ones(length(data),1); % observation probabilities
xi = zeros(size(T,1),size(T,2)); % sum of joint probabilities of adjacent transitions
if nargin > 2
    z = p;
else
    z = stat(T);
end

b(:,end) = ones(size(T,1),1);
f(:,1) = diag(sum(A{data(1)}))*z;
c(1)   = sum(f(:,1));
f(:,1) = f(:,1)/c(1);
for t = 2:length(data)
    f(:,t) = A{data(t)}*f(:,t-1);
    c(t)   = sum(f(:,t));
    f(:,t) = f(:,t)/c(t);
end

for t = length(data)-1:-1:1
    b(:,t) = b(:,t+1)'*A{data(t+1)}/c(t+1);
    xi = xi + ( (b(:,t+1)*f(:,t)') .* A{data(t+1)} )/c(t+1);
end

g = b.*f; % the gammas of the forward-backward algorithm