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
    
    
    [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
%     [err.error,resnorm,~] = ETE_integrator.step_mod(soln,err,1);
%     [err.error,resnorm,~] = ETE_integrator.step_mod(soln,err,weight);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k = 1:Error.num_iter
% %     weight = (k/Error.num_iter);
% % weight = ((Error.num_iter-k+1)/Error.num_iter);
%     weight = 0.5;%1-(1-0.5).^(Error.num_iter-k);
% 
%     % solve ETE again
%     err.error = Error.tempError;
%     ETE_integrator.em1 = 0*estError(:,err.ptr(err.M));
%     ETE_integrator.em2 = 0*estError(:,err.ptr(err.M-1));
%     ETE_integrator.pos = err.M;
%     [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
%     % modify stencil to reduce weight of new addition
% %     err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
%     err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - weight*err.error;
%     Error.R{err.count,k+1} = resnorm;
%     Error = output_error_info(Error,soln,err,initialStencil,Error.interval,err.count,k+1);
%     Error.tempError = initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1));
% end


% %% first iterative correction
% err.error = Error.tempError;
% ETE_integrator.em1 = 0*estError(:,err.ptr(err.M));
% ETE_integrator.em2 = 0*estError(:,err.ptr(err.M-1));
% ETE_integrator.pos = err.M;
% %% Normal step
% mu = err.stencil(soln.i,err.ptr(err.M+1));
% 
% [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
% err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
% s2 = var([mu,err.stencil(soln.i,err.ptr(err.M+1))],0,2);
% Error.R{err.count,2} = resnorm;
% Error = output_error_info(Error,soln,err,initialStencil,Error.interval,err.count,2);
% Error.tempError = initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1));
% weights = ones(err.N,1);
% for k = 2:Error.num_iter
%     weights(soln.i) = max(1-s2/max(s2),0.1);
%     weights(isnan(weights)) = 0.1;
%     weights(isinf(weights)) = 1;
%     % solve ETE again
%     err.error = Error.tempError;
%     ETE_integrator.em1 = 0*estError(:,err.ptr(err.M));
%     ETE_integrator.em2 = 0*estError(:,err.ptr(err.M-1));
%     ETE_integrator.pos = err.M;
%     % modify stencil to reduce weight of new addition
%     [err.error,resnorm,~] = ETE_integrator.step_mod(soln,err,weights);
%     err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
%     
%     [s2,mu] = increment_variance(s2,mu,err.stencil(soln.i,err.ptr(err.M+1)),k);
%     
%     Error.R{err.count,k+1} = resnorm;
%     Error = output_error_info(Error,soln,err,initialStencil,Error.interval,err.count,k+1);
%     Error.tempError = initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1));
% end
% 
% 
% end