d = 8;
n = 1e4;

oom = rand(3,3,d);
oom = oom./repmat(sum(sum(oom,3)),[3 1 d]);

col = 1/2*([cos(2*pi*(1:d)/d)' cos(2*pi*(1:d)/d + 2*pi/3)' cos(2*pi*(1:d)/d + 4*pi/3)'] + 1);
ent = -Inf;
while true
    oom2 = max(0,oom + 0.001*randn(3,3,d)); % small perturbation
    oom2 = oom2./repmat(sum(sum(oom2,3)),[3 1 d]);
    zs = plot_b(d,n,oom2,col);
    ent2 = nn_ent(zs');
    
    if ent2 > ent
        title(['H[B] = ' num2str(ent2)]);
        drawnow
        oom = oom2;
        ent = ent2;
    end
end
