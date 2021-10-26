function val = BDF2_res(u,um1,um2,dx,dt,nu,N,BC1,BC2)
res = -ss_residual(u,dx,nu,N);
val = u - (4/3)*um1 + (1/3)*um2 +(2/3)*dt*res;
val(1) = (u(1) - BC1);
val(N) = (u(N) - BC2);
end