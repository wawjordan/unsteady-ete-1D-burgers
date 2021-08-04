function [soln,err,ETE_integrator,Error] = advance_iterate_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,estError,Error)
% initStencil = err.stencil;
% swapStencil = err.stencil;

% [P,~] = linsolve(this.T,this.stencil(j1,this.ptr)',opts);

for k = 1:Error.num_iter
    
%     err.stencil = initStencil;
    % solve ETE again
    err.error = Error.tempError;
    ETE_integrator.em1 = 0*estError(:,err.ptr(err.M));
    ETE_integrator.em2 = 0*estError(:,err.ptr(err.M-1));
%     err.error = estError(:,err.ptr(err.M));
    ETE_integrator.pos = err.M;
    
    
    % modify stencil to reduce weight of new addition
    weight = (k/Error.num_iter).^4;
%     weight = (1-1e-2).^(Error.num_iter-k);
%     [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
%     [err.error,resnorm,~] = ETE_integrator.step_mod(soln,err,1);
    [err.error,resnorm,~] = ETE_integrator.step_mod(soln,err,weight);
    err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
    Error.R{err.count,k+1} = resnorm;
    Error = output_error_info(Error,soln,err,initialStencil,Error.interval,err.count,k+1);
%     plot(Error.tempError,'b');
    Error.tempError = initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1));
%     plot(ETE_integrator.em1,'b');
%     plot(ETE_integrator.em2,'b');
%     plot((initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1))))
% plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)))
end
% plot((initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1))),'g')
end