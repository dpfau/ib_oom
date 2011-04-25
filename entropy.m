function H = entropy(p)
% For a matrix of probability vectors, returns a column vector with the
% entropy of each one.

log_p = log(p);
log_p(isinf(log_p)) = 0; % 0log0 := 0
H = -tprod(p,[-1 2],log_p,[-1 2],'n');