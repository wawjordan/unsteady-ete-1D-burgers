function [soln,err,ETE_integrator,Error,STENCIL] = advance_ETE2(soln,err,ETE_integrator,bndry_cond,STENCIL,Error)

err.count = soln.count;
err.t(err.ptr(1)) = soln.t;

err.stencil = STENCIL(1).S;

err.ptr = circshift(err.ptr,-1);

ETE_integrator.pos = err.M;
err.error = STENCIL(1).S(soln.i,err.ptr(err.M+1)) - STENCIL(2).S(soln.i,err.ptr(err.M+1));
ETE_integrator.em1 = STENCIL(1).S(soln.i,err.ptr(err.M)) - STENCIL(2).S(soln.i,err.ptr(err.M));
ETE_integrator.em2 = STENCIL(1).S(soln.i,err.ptr(err.M-1)) - STENCIL(2).S(soln.i,err.ptr(err.M-1));

[err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
err.stencil = STENCIL(2).S;
Error.tempError = err.error;

for i = 2:length(STENCIL)
    STENCIL(i).S(soln.i,err.ptr(err.M+1)) = STENCIL(i).S(soln.i,err.ptr(err.M+1)) - err.error;
end
err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;


Error.R{err.count,1} = resnorm;
Error = output_error_info(Error,soln,err,STENCIL(1).S,Error.interval,err.count,1);

% plot((STENCIL(1).S(soln.i,err.ptr(err.M+1))-STENCIL(2).S(soln.i,err.ptr(err.M+1))),'r')

end