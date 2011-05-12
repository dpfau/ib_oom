addpath '/Users/davidpfau/Documents/Wood Group/PDIA/scripts'

seq  = loadByChar('/Users/davidpfau/Documents/Wood Group/PDIA/data/aiw.train');
seq_ = loadByChar('/Users/davidpfau/Documents/Wood Group/PDIA/data/aiw.test');

for t = 1:length(seq)
    seq{t} = seq{t} + 1;
end
for t = 1:length(seq_)
    seq_{t} = seq_{t} + 1;
end

oom = ib_oom(seq,seq_,500,5);