function [soln,err,integrator,STENCIL,Primal,Error] = advance_primal2(soln,err,integrator,STENCIL,bndry_cond,Primal,Error)
soln.count = soln.count + 1;
soln.t = soln.t + soln.dt;
fprintf('t = %f\n',soln.t);
[soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
% plot(soln.error,'k')

[Primal,Error] = output_primal_info(Primal,Error,soln,resnorm,Primal.interval,soln.count);

% update stencil
err.stencil(:,err.ptr(1)) = soln.U;
for i = 1:length(STENCIL)
    STENCIL(i).S(:,err.ptr(1)) = soln.U;
end
end