function [soln,err,ETE_integrator,estError,Error] = init_stencil_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,Error)
% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = err.M + 1;

estError = zeros(Nds,stenLength);

% solve ETE for each solution in stencil (march forward in time)
for i = 2:stenLength
    err.count = i;
    ETE_integrator.pos = i-1;
    [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
    estError(:,i) = err.error;
    Error.R{i,1} = resnorm;
end
Error.tempError = err.error;
% Correct solutions in stencil
err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
for i = 2:stenLength
Error = output_error_info(Error,soln,err,initialStencil,Error.interval,i,1);
end

end