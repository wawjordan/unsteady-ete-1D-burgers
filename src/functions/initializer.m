function [soln,integrator,R,E,t] = initializer(soln,integrator,bndry_cond,max_steps,p)
    R = cell(max_steps);
    E = nan(max_steps,soln.neq);
    t = nan(max_steps,1);
    
    soln.t=soln.t0-soln.dt;
    soln.U = soln.ExactSolution;
    soln.t=soln.t+soln.dt;
    t(1) = soln.t;
    [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
    soln.U = bndry_cond.enforce(soln,soln.U);
    soln.error = soln.U(soln.i_low:soln.i_high) - soln.ExactSolution(soln.i_low:soln.i_high);
    for j = 1:soln.neq
        E(1,j) = norm(soln.error(:,j),p)/(soln.grid.imax^(1/p));
    end
    R{1} = resnorm;
end