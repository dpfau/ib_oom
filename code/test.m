p = zeros(7,1);
s = zeros(7,1);
z = zeros(length(seq),1);
z(1) = 7;

d = ceil(rand(2,7)*7);
for i = 2:length(seq)
    z(i) = d(seq(i-1),z(i-1));
end
for i = 1:7
    s(i) = sum(z == i)/length(seq);
    if s(i) == 0
        p(i) = 1;
    else
        p(i) = sum(z == i & seq == 1)/sum(z == i);
    end
end
oom = zeros(7,7,2);
for i = 1:7
    oom(d(1,i),i,1) = p(i);
    oom(d(2,i),i,2) = 1 - p(i);
end
[a,b,c,l,xi] = for_back_oom(seq,oom);
x = entropy(s) - mean(entropy(normalize(diag(s)*b)))
y = entropy(s) - mean(entropy(a)) - mean(entropy(normalize(diag(s)*b))) + mean(entropy(c))
z = sum(log((seq == 1) .* p(z) + (seq == 2) .* (1 - p(z))))/length(seq)