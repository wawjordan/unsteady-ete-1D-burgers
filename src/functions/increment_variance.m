function [s2,mu] = increment_variance(s2,mu,x,N)

mu_old = mu;
mu = (1/N)*(x+(N-1)*mu_old);
s2_old = s2;
s2 = ((N-1)*s2_old + N*(mu_old - mu).^2 + (x - mu).^2)/N;

end