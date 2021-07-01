function Error = output_error_info(Error,soln,err,initialStencil,out_interval,i,iter)
p = [1,2,inf];
N = length(p);
M = err.M;
if i <= M+1
    current_error = initialStencil(soln.i,i) ...
        - soln.calc_exact(soln.grid.x(soln.i),err.t(i));
    estimated_error = initialStencil(soln.i,i)-err.stencil(soln.i,i);
else
% if i <= M+1
%     current_error = initialStencil(soln.i,err.ptr(i)) ...
%         - soln.calc_exact(soln.grid.x(soln.i),err.t(err.ptr(i)));
%     estimated_error = initialStencil(soln.i,err.ptr(i))-err.stencil(soln.i,err.ptr(i));
% else
current_error = initialStencil(soln.i,err.ptr(M+1)) ...
        - soln.calc_exact(soln.grid.x(soln.i),err.t(err.ptr(M+1)));
    estimated_error = initialStencil(soln.i,err.ptr(M+1))-err.stencil(soln.i,err.ptr(M+1));
end
    for k = 1:N
        Error.E(i,iter,k) = norm(estimated_error - current_error,p(k))/(soln.grid.imax^(1/p(k)));
    end
    Error.Et(:,iter,1) = Error.Et(:,iter,1) + abs(estimated_error - current_error);
    Error.Et(:,iter,2) = Error.Et(:,iter,2) + (estimated_error - current_error).^2;
    Error.Et(:,iter,3) = max(Error.Et(:,iter,3),abs(estimated_error - current_error));
    
    if mod(i,out_interval)==0
        Error.out.Eerror{i,iter} = estimated_error-current_error;
        Error.out.error{i,iter} = estimated_error;
%         Error.out.Eerror{i} = estimated_error-current_error;
%         Error.out.error{i} = estimated_error;
    end
end