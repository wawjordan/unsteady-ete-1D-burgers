function [soln,err,ETE_integrator,Error] = advance_iterate_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,estError,Error)
% initStencil = err.stencil;
% swapStencil = err.stencil;
for k = 1:Error.num_iter
    
%     err.stencil = initStencil;
    % solve ETE again
    err.error = Error.tempError;
    ETE_integrator.em1 = 0*estError(:,err.ptr(err.M));
    ETE_integrator.em2 = 0*estError(:,err.ptr(err.M-1));
%     err.error = estError(:,err.ptr(err.M));
    ETE_integrator.pos = err.M;
    [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
    err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
    Error.R{err.count,k+1} = resnorm;
    Error = output_error_info(Error,soln,err,initialStencil,Error.interval,err.count,k+1);
    Error.tempError = initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1));
% plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)))
end
% plot((initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1))))
end