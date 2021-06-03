function [soln,integrator,Primal] = solver2(soln,integrator,bndry_cond,max_steps)
    p = 1;
    M = 2;
    out_interval = max_steps;
    
    Primal.R = cell(max_steps+M,soln.neq);
    Primal.E = nan(max_steps+M,soln.neq);
    Primal.Et = zeros(soln.grid.imax,soln.neq);
    Primal.Etf = 0;
    Primal.t = nan(max_steps+M,1);
    Primal.out.t = nan(max_steps+M,1);
    Primal.out.error = cell(max_steps+M,soln.neq);
    Primal.out.u = cell(max_steps+M,soln.neq);
    
    %% Initial Conditions
    soln.t=soln.t0-(M+1)*soln.dt;
    for i = 1:M
        soln.t=soln.t+soln.dt;
        soln.count = i;
        soln.U = soln.ExactSolution;
        for j = 1:soln.neq
            Primal.E(i,j) = 0;
        end
        Primal.R{i} = 0;
        if mod(i,out_interval)==0
            Primal.out.t(i) = soln.t;
            Primal.out.error{i} = soln.U(soln.i)-soln.ExactSolution(soln.i);
            Primal.out.u{i} = soln.U(soln.i);
        end
    end
    while (soln.count < max_steps+M)&& (soln.t <= soln.tf)
        i = i+1;
        soln.t=soln.t+soln.dt;
        soln.count = i;
        Primal.t(i) = soln.t;
        [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
        soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);

        for j = 1:soln.neq
            Primal.E(i,j) = norm(soln.error(:,j),p)/(soln.grid.imax^(1/p));
            if p == 1
                Primal.Et(:,j) = Primal.Et(:,j) + abs(soln.error(:,j));
            elseif p == 2
                Primal.Et(:,j) = Primal.Et(:,j) + soln.error(:,j).^2;
            elseif isinf(p)
                Primal.Et(:,j) = max(Primal.Et(:,j),abs(soln.error(:,j)));
            end
        end
        Primal.R{i} = resnorm;
        fprintf('.');
        if mod(i,10)==0
            fprintf('\n');
        end
        if mod(i,out_interval)==0
            Primal.out.t(i) = soln.t;
            Primal.out.error{i} = soln.U(soln.i)-soln.ExactSolution(soln.i);
            Primal.out.u{i} = soln.U(soln.i);
        end
    end
    fprintf('\n');
    Primal.out.t(i) = soln.t;
    Primal.out.error{i} = soln.U(soln.i)-soln.ExactSolution(soln.i);
    Primal.out.u{i} = soln.U(soln.i);

    Primal.R = Primal.R(~cellfun('isempty',Primal.R));
    Primal.R = cell2mat(Primal.R);
    Primal.E = Primal.E(~isnan(Primal.t),:);
    Primal.Et = Primal.Et/((soln.count-M)^(1/p));
    for j = 1:soln.neq
        Primal.Etf(1,j) = norm(Primal.Et(:,j),p)/(soln.grid.imax^(1/p));
    end
    Primal.t = Primal.t(~isnan(Primal.t));
    Primal.out.t = Primal.out.t(~isnan(Primal.out.t));
    Primal.out.error = Primal.out.error(~cellfun('isempty',Primal.out.error));
    Primal.out.u = Primal.out.u(~cellfun('isempty',Primal.out.u));
end