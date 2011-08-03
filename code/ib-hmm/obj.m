function x = obj(gamma,gam_tot,xi,O,T,beta)

a = log(O(:));
a(isinf(a)) = 0;

b = log(stat(join_oom(O,T)));
b(isinf(b)) = 0;

c = log(T(:));
c(isinf(c)) = 0;

x = beta * gamma(:)' * a - gam_tot' * b + xi(:)' * c;