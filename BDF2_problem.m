function val = BDF2_problem(u,um1,um2,dx,dt,nu,N,BC1,BC2)
%     u(1) = BC1;
%     u(N) = BC2;
    val = BDF2_LHS(u,um1,um2) - BDF2_RHS(u,dx,dt,nu,N);
end