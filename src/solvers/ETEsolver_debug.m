function [soln,err,integrator,ETE_integrator,Primal,Error] = ETEsolver_debug(soln,err,integrator,ETE_integrator,bndry_cond,max_steps)
    p = 2;
    M = err.M;
    
    Primal.R = cell(max_steps+M,soln.neq);
    Primal.E = nan(max_steps+M,soln.neq);
    Primal.t = nan(max_steps+M,1);
    Primal.exact_error = cell(max_steps+M,1);
    Primal.tau_exact = cell(max_steps+M,1);
    Primal.u_exact = cell(max_steps+M,1);
    Primal.u = cell(max_steps+M,1);

    Error.R = cell(max_steps,soln.neq);
    Error.E = nan(max_steps,soln.neq);
    Error.t = nan(max_steps,1);
    Error.est_error = cell(max_steps,1);
    Error.stencil = cell(max_steps,1);
    Error.tau_est = cell(max_steps,1);
    
    %% Initial Conditions
    
    soln.t=soln.t0-(M/2+1)*soln.dt;
    soln.U = soln.ExactSolution;
    temp_sol = soln.U;
    for i = 1:M
        soln.t=soln.t+soln.dt;
        fprintf('soln.t = %0.4f\n',soln.t)
        Primal.t(i) = soln.t;
        [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
        if i > 1
            Primal.tau_exact{i} = soln.residual(soln.ExactSolution) + (soln.U(soln.i)-err.stencil(soln.i,i-1))/soln.dt;
        else
            Primal.tau_exact{1} = soln.residual(soln.ExactSolution) + (soln.U(soln.i)-temp_sol(soln.i))/soln.dt;
        end
        
        soln.count = i;
        soln.error = soln.U(soln.i_low:soln.i_high) - soln.ExactSolution(soln.i_low:soln.i_high);
        Primal.exact_error{i} = soln.error;
        Primal.u_exact{i} = soln.ExactSolution;
        Primal.u{i} = soln.U;
        for j = 1:soln.neq
            Primal.E(i,j) = norm(soln.error(:,j),p)/(soln.grid.imax^(1/p));
        end
        Primal.R{i} = resnorm;
        err = err.update_time_stencil(soln);
    end
    err.count = 0;
    while (err.count < max_steps)&&(err.t(err.ptr(M/2+1)) <= err.tf)
        i = i+1;
        soln.t=soln.t+soln.dt;
        fprintf('soln.t = %0.4f\n',soln.t)
        Primal.t(i) = soln.t;
        [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
        soln.count = i;
        soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
        
        for j = 1:soln.neq
            Primal.E(i,j) = norm(soln.error(:,j),p)/(soln.grid.imax^(1/p));
        end
        Primal.R{i} = resnorm;
        err = err.update_time_stencil(soln);
        err.count = err.count+1;
        fprintf('err.t = %0.4f\n',err.t(err.ptr(M/2+1)))
        Primal.exact_error{i} = soln.error;
        Primal.u_exact{i} = soln.ExactSolution;
        Primal.u{i} = soln.U;
        Primal.tau_exact{i} = soln.residual(soln.ExactSolution) + (soln.U(soln.i)-Primal.u{i-1}(soln.i))/soln.dt;

        Error.t(err.count) = err.t(err.ptr(M/2+1));
        current_error = err.stencil(soln.i,err.ptr(1)) - soln.calc_exact(soln.grid.x(soln.i),err.t(err.ptr(1)));
        [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
        Error.E(err.count,1) = norm(err.error - current_error,p)/(soln.grid.imax^(1/p));
        Error.R{err.count} = resnorm;
        Error.est_error{err.count} = err.error;
        Error.stencil{err.count} = err.stencil(:,err.ptr);
        Error.tau_est{err.count} = ETE_integrator.tau;
        
    end
    
    Primal.R = Primal.R(~cellfun('isempty',Primal.R));
    Primal.R = cell2mat(Primal.R);
    Primal.E = Primal.E(~isnan(Primal.t),:);
    Primal.t = Primal.t(~isnan(Primal.t));
    
    Primal.exact_error = Primal.exact_error(~cellfun('isempty',Primal.exact_error));
    Primal.tau_exact = Primal.tau_exact(~cellfun('isempty',Primal.tau_exact));
    Primal.u_exact = Primal.u_exact(~cellfun('isempty',Primal.u_exact));
    Primal.u = Primal.u(~cellfun('isempty',Primal.u));
    
    Error.R = Error.R(~cellfun('isempty',Error.R));
    Error.R = cell2mat(Error.R);
    
    Error.est_error = Error.est_error(~cellfun('isempty',Error.est_error));
    Error.stencil = Error.stencil(~cellfun('isempty',Error.stencil));
    Error.tau_est = Error.tau_est(~cellfun('isempty',Error.tau_est));
    
    Error.E = Error.E(~isnan(Error.t),:);
    Error.t = Error.t(~isnan(Error.t));
end