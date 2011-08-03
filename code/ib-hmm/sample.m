function x = sample(p)
% return random sample from probability vector p

cusum = cumsum(p);
x = find( cusum > rand*cusum(end), 1 );