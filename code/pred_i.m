function x = pred_i(f,b,g,z)
% Gives an estimate of I[ X_{-\infty:0}; X_{1:\infty} ] = 
% I[ X_{-\infty:0}; Z_1 ] - I[ X_{-\infty:0}; Z_1 | X_{1:\infty} ] = 
% H[ Z_1 ] - H[ Z_1 | X_{-\infty:0} ] - H[ Z_1 | X_{1:\infty} ] + H[ Z_1 | X_{-\infty:\infty} ]

b2 = b .* (z*ones(1,size(b,2)));
b2 = (ones(size(b,1),1)*(1./sum(b2))) .* b2;

x = entropy(z) - mean(entropy(f)) - mean(entropy(b2)) + mean(entropy(g));