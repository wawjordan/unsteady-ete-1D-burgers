classdef back_diff_2mod2 < time_integrator_type
    properties
        imax, i_low, i_high, i
        Rnorm, Rinit
        um1, um2
        start
    end
    methods
        function this = back_diff_2mod2(soln)
            this.i_low = soln.i_low;
            this.i_high = soln.i_high;
            this.imax = soln.grid.imax;
            this.i = soln.i;
            this.um1 = soln.calc_exact(soln.grid.x,soln.t - min(soln.dt));
            this.um2 = soln.calc_exact(soln.grid.x,soln.t - 2*min(soln.dt));
            this.start = true;
        end
        function [u_new,R,this] = step(this,soln,bndry_cond)
            u_new = soln.U;
            uex = soln.ExactSolution(this.i);
            u_new = bndry_cond.enforce(soln,u_new);
            F  = @(u) BDF2_res(u,this.um1(this.i),this.um2(this.i),soln.grid.dx(this.i),soln.dt,soln.nu,soln.grid.imax,uex(1),uex(end));
            dF = @(u) BDF2_jac(u,soln.grid.dx(this.i),soln.dt,soln.nu,soln.grid.imax);
%             [xk,Rks] = newton_with_backtracking(x0,F,dF,gamma,res_tol,max_iter)
            [u_new(this.i),R] = newton_with_backtracking(u_new(this.i),F,dF,0.5,this.newton_tol,this.newton_max_iter);
            this.um2 = this.um1;

            
            this.um1 = u_new;
            
        end
    end
end