function seq = gen(len,O,T)

if nargin == 2
    [O,T] = split_oom(O);
end
seq = zeros(len,1);
z = sample(stat(join_oom(O,T)));
for i = 1:len
    seq(i) = sample(O(:,z));
    z = sample(T(:,z,seq(i)));
end