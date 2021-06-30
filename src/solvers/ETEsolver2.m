function [soln,err,integrator,ETE_integrator,Primal,Error] = ETEsolver2(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval)
    err = err.reset;
    p = [1,2,inf];
    N = length(p);
%     out_interval = max_steps;
    
    % reconstruction stencil size
    M = err.M;
    
    Primal.R = cell(max_steps+M,soln.neq);
    Primal.E = nan(max_steps+M,soln.neq,N);
    Primal.Ef = zeros(1,soln.neq,N);
    Primal.Et = zeros(soln.grid.imax,soln.neq,N);
    Primal.Etf = zeros(1,soln.neq,N);
    Primal.t = nan(max_steps+M,1);
    Primal.out.t = nan(max_steps+M,1);
    Primal.out.error = cell(max_steps+M,soln.neq);
    Primal.out.u = cell(max_steps+M,soln.neq);
    
    Error.R = cell(max_steps,soln.neq);
    Error.E = nan(max_steps,soln.neq,N);
    Error.Ef = zeros(1,soln.neq,N);
    Error.Et = zeros(soln.grid.imax,soln.neq,N);
    Error.Etf = zeros(1,soln.neq,N);
    Error.t = nan(max_steps,1);
    Error.out.t = nan(max_steps,1);
    Error.out.Eerror = cell(max_steps,soln.neq);
    Error.out.error = cell(max_steps,soln.neq);
    
    %% Initial Conditions
%   soln.t=soln.t0-(M/2+1)*soln.dt;
%   changing to 1-sided stencil
    soln.t=soln.t0-(M+1)*soln.dt;
    for i = 1:M
        soln.t=soln.t+soln.dt;
        soln.count = i;
        soln.U = soln.ExactSolution;
        for j = 1:soln.neq
            Primal.E(i,j) = 0;
        end
        Primal.R{i} = 0;
        err = err.update_time_stencil(soln);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    err.count = 0;
    while (err.count < max_steps)&&(err.t(err.ptr(M)) < err.tf)
        i = i+1;
        soln.t=soln.t+soln.dt;
        soln.count = i;
        Primal.t(i) = soln.t;
        [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
        soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);

        for j = 1:soln.neq
            for k = 1:N
                Primal.E(i,j,k) = norm(soln.error(:,j),p(k))/(soln.grid.imax^(1/p(k)));
            end
        end
        Primal.Et(:,j,1) = Primal.Et(:,j,1) + abs(soln.error(:,j));
        Primal.Et(:,j,2) = Primal.Et(:,j,2) + soln.error(:,j).^2;
        Primal.Et(:,j,3) = max(Primal.Et(:,j,3),abs(soln.error(:,j)));
        Primal.R{i} = resnorm;
        if mod(i-M-1,out_interval)==0
            Primal.out.t(i) = soln.t;
            Primal.out.error{i} = soln.U(soln.i)-soln.ExactSolution(soln.i);
            Primal.out.u{i} = soln.U(soln.i);
        end
        % Here we have a saved copy of the uh0 solution
        % now we can operate on the soln.U array and recover the original
        % primal solution later
        Uh0 = soln.U;
        
%=========================================================================%
        err = err.update_time_stencil(soln);
        err.count = err.count+1;
        fprintf('Iter - %0.4d\n',err.count);
        Error.t(err.count) = err.t(err.ptr(M+1));
        
        Ncorr = 1;
        for nc = 1:Ncorr
            [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
%             err.stencil(soln.i,err.ptr(M+1)) = err.stencil(soln.i,err.ptr(M+1)) - err.error;
            soln.U(soln.i) = soln.U(soln.i) - err.error;
        end
        err.error = Uh0(soln.i) - soln.U(soln.i);

%         err.error = Uh0(soln.i) - err.stencil(soln.i,err.ptr(M+1));
%         err.stencil(:,err.ptr(M+1)) = Uh0;




        soln.U = Uh0;
%         for nc = 1:Ncorr
%             [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
%             err.stencil(soln.i,err.ptr(M+1)) = err.stencil(soln.i,err.ptr(M+1)) - err.error;
%         end
%         err.error = Uh0(soln.i) - err.stencil(soln.i,err.ptr(M+1));
%=========================================================================%
        current_error = err.stencil(soln.i,err.ptr(M+1)) ...
            - soln.calc_exact(soln.grid.x(soln.i),err.t(err.ptr(M+1)));
        
        for j = 1:soln.neq
            for k = 1:N
                Error.E(err.count,j,k) = norm(err.error(:,j) - current_error(:,j),p(k))/(soln.grid.imax^(1/p(k)));
            end
        end
        Error.Et(:,j,1) = Error.Et(:,j,1) + abs(err.error(:,j) - current_error(:,j));
        Error.Et(:,j,2) = Error.Et(:,j,2) + (err.error(:,j) - current_error(:,j)).^2;
        Error.Et(:,j,3) = max(Error.Et(:,j,3),abs(err.error(:,j) - current_error(:,j)));
        Error.R{err.count} = resnorm;
        if mod(err.count-1,out_interval)==0&&(Error.t(err.count)-1e-6<=soln.tf)
            Error.out.t(err.count) = Error.t(err.count);
            Error.out.Eerror{err.count} = err.error-current_error;
            Error.out.error{err.count} = err.error;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indP = find(abs(Primal.t-soln.tf)<1e-6,1);
    % Primal.R = Primal.R(~cellfun('isempty',Primal.R));
    Primal.R = Primal.R(M+1:indP);
    Primal.R = cell2mat(Primal.R);
    % Primal.t = Primal.t(~isnan(Primal.t));
    Primal.t = Primal.t(M+1:indP);
    % Primal.E = Primal.E(~isnan(Primal.t),:,:);
    Primal.E = Primal.E(M+1:indP,:,:);
    for k = 1:N
        if ~isinf(p(k))
            Primal.Et(:,:,k) = (Primal.Et(:,:,k)/(indP-M-1)).^(1/p(k));
        end
    end
    
    for j = 1:soln.neq
        for k = 1:N
            Primal.Etf(1,j,k) = norm(Primal.Et(:,j,k),p(k))/(soln.grid.imax^(1/p(k)));
            %Primal.Ef(1,j,k) = norm(Primal.E(M/2:indP,j,k),p(k))/((indP-M/2)^(1/p(k)));
            Primal.Ef(1,j,k) = norm(Primal.E(:,j,k),p(k))/((indP-M)^(1/p(k)));
        end
    end
    Primal.out.t = Primal.out.t(~isnan(Primal.out.t));
    Primal.out.error = Primal.out.error(~cellfun('isempty',Primal.out.error));
    Primal.out.u = Primal.out.u(~cellfun('isempty',Primal.out.u));
    
    
    indE = find(abs(Error.t-soln.tf)<1e-6,1);
    % Error.R = Error.R(~cellfun('isempty',Error.R));
    Error.R = Error.R(1:indE);
    Error.R = cell2mat(Error.R);
    % Error.t = Error.t(~isnan(Error.t));
    Error.t = Error.t(1:indE);
    % Error.E = Error.E(~isnan(Error.t),:,:);
    Error.E = Error.E(1:indE,:,:);
    for k = 1:N
        if ~isinf(p(k))
            Error.Et(:,:,k) = (Error.Et(:,:,k)/err.count).^(1/p(k));
        end
    end
    for j = 1:soln.neq
        %Error.Ef(1,j,:) = Error.E(ind,j,:);
        for k = 1:N
            Error.Etf(1,j,k) = norm(Error.Et(:,j,k),p(k))/(soln.grid.imax^(1/p(k)));
            Error.Ef(1,j,k) = norm(Error.E(1:indE,j,k),p(k))/(indE^(1/p(k)));
        end
    end
    Error.out.t = Error.out.t(~isnan(Error.out.t));
    Error.out.Eerror = Error.out.Eerror(~cellfun('isempty',Error.out.Eerror));
    Error.out.error = Error.out.error(~cellfun('isempty',Error.out.error));
    
    
    Primal.out.t = Primal.out.t(1:length(Error.out.t));
    Primal.out.error = Primal.out.error(1:length(Error.out.t));
    Primal.out.u = Primal.out.u(1:length(Error.out.t));
end