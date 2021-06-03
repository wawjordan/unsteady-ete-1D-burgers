classdef exact_BC < boundary_condition
    methods
        function  u = enforce(this,soln,u)
            u(1:this.i_low-1) = soln.calc_exact(soln.grid.x(1:this.i_low-1),soln.t);
            u(this.i_high+1:end) = soln.calc_exact(soln.grid.x(this.i_high+1:end),soln.t);
        end
    end
end