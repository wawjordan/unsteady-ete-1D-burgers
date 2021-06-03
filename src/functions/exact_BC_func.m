function u = exact_BC(soln,u)
u(1:soln.i_low-1) = soln.ExactSolution(1:soln.i_low-1);
    u(soln.i_high+1:end) = soln.ExactSolution(soln.i_high+1:end);
end

