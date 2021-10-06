function u = extrap_BC(soln,u)
for i = soln.i_low-1:-1:1
    u(i) = 2*u(i+1) - u(i+2);
end
for i = soln.i_high+1:length(u)
    u(i) = 2*u(i-1) - u(i-2);
end