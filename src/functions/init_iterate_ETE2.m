function [soln,err,ETE_integrator,Error,STENCIL] = init_iterate_ETE2(soln,err,ETE_integrator,bndry_cond,STENCIL,Error)
% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = err.M + 1;

estError = zeros(Nds,stenLength);

for k = 1:Error.num_iter
    % solve ETE again
    err.error = 0*err.error;
    if (isa(ETE_integrator,'back_diff_2_ETE')||isa(ETE_integrator,'back_diff_2_ETEmod')||isa(ETE_integrator,'back_diff_2_ETEmod2'))
        ETE_integrator.em1 = 0*err.error;
        ETE_integrator.em2 = 0*err.error;
    end
    for i = 2:stenLength
        err.count = i;
        ETE_integrator.pos = i-1;
        [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
        estError(:,i) = err.error;
        Error.R{i,k+1} = resnorm;
    end
    
    % Correct solutions in stencil
    err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
    STENCIL(k+2).S = err.stencil;
    for i = 2:stenLength
        Error = output_error_info(Error,soln,err,STENCIL(1).S,Error.interval,i,k+1);
    end
%     plot((STENCIL(1).S(soln.i,:)-err.stencil(soln.i,:)))
end

end