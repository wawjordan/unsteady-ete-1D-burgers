function RHS = BDF2_RHS(u,dx,dt,nu,N)
    RHS = (2/3)*dt*ss_residual(u,dx,nu,N);
end