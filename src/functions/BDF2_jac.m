function val = BDF2_jac(u,dx,dt,nu,N)
val = zeros(N,3);
for ii = 2:N-1
    val(ii,1) = (nu/dx(ii-1)^2) + 0.5*u(ii-1)/dx(ii-1);
    val(ii,2) = -2*(nu/dx(ii)^2);
    val(ii,3) = (nu/dx(ii+1)^2) - 0.5*u(ii+1)/dx(ii+1);
end
val = -(2/3)*dt*val;  % modify jacobian for BDF2 algorithm
val(:,2) = val(:,2) + 1;
val(1,2) = 1;
val(end,2) = 1;
end