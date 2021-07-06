function [soln,err,integrator,initialStencil,Primal,Error] = advance_primal(soln,err,integrator,initialStencil,bndry_cond,Primal,Error)
soln.count = soln.count + 1;
soln.t = soln.t + soln.dt;
[soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
plot(soln.error,'k')

[Primal,Error] = output_primal_info(Primal,Error,soln,resnorm,Primal.interval,soln.count);

% update stencil
err.stencil(:,err.ptr(1)) = soln.U;
initialStencil(:,err.ptr(1)) = soln.U;
end