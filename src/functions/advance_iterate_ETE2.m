function [soln,err,ETE_integrator,Error,STENCIL] = advance_iterate_ETE2(soln,err,ETE_integrator,bndry_cond,STENCIL,Error)
for k = 1:Error.num_iter
    % solve ETE again
    % use error estimates using current iteration level
    err.error = STENCIL(1).S(soln.i,err.ptr(err.M+1)) - STENCIL(k+2).S(soln.i,err.ptr(err.M+1));
    if (isa(ETE_integrator,'back_diff_2_ETE')||isa(ETE_integrator,'back_diff_2_ETEmod')||isa(ETE_integrator,'back_diff_2_ETEmod2'))
        ETE_integrator.em1 = STENCIL(1).S(soln.i,err.ptr(err.M)) - STENCIL(k+2).S(soln.i,err.ptr(err.M));
        ETE_integrator.em2 = STENCIL(1).S(soln.i,err.ptr(err.M-1)) - STENCIL(k+2).S(soln.i,err.ptr(err.M-1));
    end
    ETE_integrator.pos = err.M;
    
    % reset stencil to uncorrected stencil (prevents correction from being
    % applied twice)
    err.stencil = STENCIL(1).S;
    % solve ETE
    [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
    
    % correct stencil
    STENCIL(k+2).S(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
    
    % update err.stencil for use with "output_error_info" function
    err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
    
    % save iteration info
    Error.R{err.count,k+1} = resnorm;
    Error = output_error_info(Error,soln,err,STENCIL(1).S,Error.interval,err.count,k+1);
    
    % debug/local error plots
%     plot(err.error,'b');
%     plot(ETE_integrator.em1,'b');
%     plot(ETE_integrator.em2,'b');
%     plot((STENCIL(1).S(soln.i,err.ptr(err.M+1))-STENCIL(k+2).S(soln.i,err.ptr(err.M+1))))
end
% plot((STENCIL(1).S(soln.i,err.ptr(err.M+1))-STENCIL(k+2).S(soln.i,err.ptr(err.M+1))),'g')

end