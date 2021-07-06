function [soln,err,integrator,exactError,initialStencil,Primal,Error] = fill_stencil(soln,err,integrator,bndry_cond,Primal,Error)
% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = err.M + 1;
exactError = zeros(Nds,stenLength);
for i = 2:stenLength
    soln.count = i;
    soln.t = soln.t + soln.dt; % (update time for correct application of exact BCs)
    [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
    err.stencil(:,i) = soln.U;
    err.t(i) = soln.t;
    soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
    exactError(:,i) = soln.error;
    [Primal,Error] = output_primal_info(Primal,Error,soln,resnorm,Primal.interval,i);
end
initialStencil = err.stencil;
end