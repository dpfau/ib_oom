function [ent zs] = ent_b(oom)
% for a given oom (or oom representation of an hmm), approximate H[B] by
% generating data and using the nearest neighbors entropy estimate

zs = zeros(size(oom,1),1e5);

z = stat(oom);
p_xz = squeeze(sum(oom))';
seq = zeros(1e5,1);
for i = 1:1e5
    prob = p_xz*z;
    cusum = cumsum(prob);
    seq(i) = find(cusum >= cusum(end)*rand,1);
    z = oom(:,:,seq(i))*z;
    z = z/sum(z);
    zs(:,i) = z;
end

ent = nn_ent(zs');

x = [0.8165, -0.4082, -0.4082; 0, 0.7071, -0.7071];
proj = x*(zs-1/3*ones(size(zs))); % project from the 3-simplex to a triangle

clf
line([x(2,1) x(2,1) x(2,2); x(2,2) x(2,3) x(2,3)],[x(1,1) x(1,1) x(1,2); x(1,2) x(1,3) x(1,3)],'Color','k'); 
axis square; hold on

for i = 1:size(oom,3)
    scatter(proj(2,seq==i),proj(1,seq==i),1,rand(1,3));
end
title(['H[B] = ' num2str(ent)]);