function [soln,err,ETE_integrator,estError,Error] = advance_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,estError,Error)

err.count = soln.count;
swapStencil = err.stencil;
err.stencil = initialStencil;
err.t(err.ptr(1)) = soln.t;
err.ptr = circshift(err.ptr,-1);

ETE_integrator.pos = err.M;
% ETE_integrator.em1 = Error.tempError;
if (isa(ETE_integrator,'back_diff_2_ETE')||isa(ETE_integrator,'back_diff_2_ETEmod')||isa(ETE_integrator,'back_diff_2_ETEmod2'))
    ETE_integrator.em1 = estError(:,err.ptr(err.M));
    ETE_integrator.em2 = estError(:,err.ptr(err.M-1));
end
% err.error = Error.tempError;
err.error = estError(:,err.ptr(err.M));
[err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);

Error.tempError = err.error;
err.stencil = swapStencil;
err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
Error.R{err.count,1} = resnorm;
Error = output_error_info(Error,soln,err,initialStencil,Error.interval,err.count,1);
% plot(err.error,'r')

estError(:,err.ptr(err.M+1))= err.error;

% plot((initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1))),'r')
end