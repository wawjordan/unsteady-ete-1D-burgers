function [soln,err,ETE_integrator,Error] = init_iterate_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,Error)
% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = err.M + 1;

estError = zeros(Nds,stenLength);

for k = 1:Error.num_iter
    % solve ETE again
    err.error = 0*err.error;
    ETE_integrator.em1 = 0*err.error;
    ETE_integrator.em2 = 0*err.error;
    for i = 2:stenLength
        err.count = i;
        ETE_integrator.pos = i-1;
        [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
        estError(:,i) = err.error;
        Error.R{i,k+1} = resnorm;
    end
    
    % Correct solutions in stencil
    err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
    for i = 2:stenLength
        Error = output_error_info(Error,soln,err,initialStencil,Error.interval,i,k+1);
    end
    
    plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)))
end

end