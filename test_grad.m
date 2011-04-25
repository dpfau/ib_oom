oom = betarnd(.1, 10, 10, 10, 5);
oom = oom./repmat(sum(sum(oom,3)), [10 1 5]);

z_0 = zeros(10,1);
z_0(1) = 1;

seq = generate(oom,z_0,2);

oom2 = betarnd(.1, 10, 10, 10, 5);
oom2 = oom2./repmat(sum(sum(oom2,3)), [10 1 5]);

b = 1;

dObj = zeros(10,10,5); dObj2 = zeros(10,10,5); dx = zeros(10,10,5);

for i = 1:10
    for j = 1:10
        for k = 1:5
            dx(i,j,k) = 1e-6;
            dObj(i,j,k)  = 1e6*(objective(seq, oom + dx, z_0, b) - objective(seq, oom, z_0, b));
            dObj2(i,j,k) = 1e6*(objective(seq, oom2 + dx,z_0, b) - objective(seq, oom2,z_0, b));
            dx(i,j,k) = 0;
        end
    end
end

grad  = gradient(seq, oom,  z_0, b);
grad2 = gradient(seq, oom2, z_0, b);

