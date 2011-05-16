function [obj ll ent] = score(seq,oom)

if strcmp(class(seq),'double')
    seq = {seq};
end

len = 0;
for i = 1:length(seq)
    len = len + length(seq{i});
end

zs = zeros(size(oom,1),len);

t = 0;
ll = 0; % Log likelihood of the data
for i = 1:length(seq)
    z = stat(oom); % initialize state vector
    for j = 1:length(seq{i})
        t = t+1;
        z = oom(:,:,seq{i}(j))*z;
        p = sum(z);
        ll = ll + log(p);
        z = z/p;
        zs(:,t) = z;
    end
end

ent = nn_ent(zs'); % empirical nearest-neighbors estimate of entropy of the 
                   % stationary distribution over belief vectors
obj = ent - ll;    % minimize this!