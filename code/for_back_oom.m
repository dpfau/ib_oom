 function [f b g c xi] = for_back_oom(data,oom,p)
% Given data, observation matrix O and transition matrix T, returns the
% forward vectors f, backward vectors b, and the marginal vectors g (gamma)
% as well as the probability of observations c

f = zeros(size(oom,1),length(data)); % forward probabilities
b = zeros(size(oom,1),length(data)); % backward probabilities
c = ones(length(data),1); % observation probabilities
xi = zeros(size(oom)); % sum of joint probabilities of adjacent transitions
if nargin > 2
    z = p;
else
    z = stat(oom);
end

b(:,end) = ones(size(oom,1),1);
f(:,1) = diag(sum(oom(:,:,data(1))))*z;
c(1)   = sum(f(:,1));
f(:,1) = f(:,1)/c(1);
for t = 2:length(data)
    f(:,t) = oom(:,:,data(t))*f(:,t-1);
    c(t)   = sum(f(:,t));
    f(:,t) = f(:,t)/c(t);
end

for t = length(data)-1:-1:1
    b(:,t) = b(:,t+1)'*oom(:,:,data(t+1))/c(t+1);
    xi(:,:,data(t+1)) = xi(:,:,data(t+1)) + ( (b(:,t+1)*f(:,t)') .* oom(:,:,data(t+1)) )/c(t+1);
end

g = b.*f; % the gammas of the forward-backward algorithm