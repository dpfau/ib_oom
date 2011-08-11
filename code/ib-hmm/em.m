function [o,t] = em(dat,n,b)

k = max(cellfun(@max,dat));

o = rand(n,k);
o = o./(sum(o,2)*ones(1,k));

t = rand(n,n,k);
t = t./repmat(sum(t),[n 1 1]);

while 1
    %% E Step
    g = zeros(size(o));
    gt = zeros(n,1);
    xi = zeros(size(t));
    for i = 1:length(dat)
        [a,b,c] = fwdback(dat{i},o,t);
        [g0, gt0, xi0] = suff_stat(dat{i},a,b,c,o,t);
        g = g + g0;
        gt = gt + gt0;
        xi = xi + xi0;
    end
    %% M Step
    [o,t,~] = m_step(g,gt,xi,b,o,t,stat(join_oom(o,t)));
    fprintf('Objective: %d\n',obj(g,gt,xi,o,t,b));
end