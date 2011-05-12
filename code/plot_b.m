function [zs seq oom] = plot_b(dim,len)

oom = rand(3,3,dim);
oom = oom./repmat(sum(sum(oom,3)),[3 1 dim]);

seq = generate(oom,len);
zs = zeros(3,len);
zs(:,1) = oom(:,:,seq(1))*stat(oom);
zs(:,1) = zs(:,1)/sum(zs(:,1));

for i = 2:len
    zs(:,i) = oom(:,:,seq(i))*zs(:,i-1);
    zs(:,i) = zs(:,i)/sum(zs(:,i));
end

x = [0.8165, -0.4082, -0.4082; 0, 0.7071, -0.7071];
proj = x*(zs-1/3*ones(size(zs))); % project from the 3-simplex to a triangle

clf
line([x(2,1) x(2,1) x(2,2); x(2,2) x(2,3) x(2,3)],[x(1,1) x(1,1) x(1,2); x(1,2) x(1,3) x(1,3)],'Color','k'); 
axis square; hold on

for i = 1:dim
    scatter(proj(2,seq==i),proj(1,seq==i),1,rand(1,3));
end

title(['Sensitivity: ' num2str(log(cover(zs,1e-3))/-log(1e-3))])