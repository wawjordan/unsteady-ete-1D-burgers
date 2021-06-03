classdef explicit_euler < time_integrator_type
    properties
        imax
        i_low
        i_high
        i
        nu
        dx
    end
    methods
        function this = explicit_euler(soln)
            this.i_low = soln.i_low;
            this.i_high = soln.i_high;
            this.imax = soln.grid.imax;
            this.i = soln.i;
            this.nu = soln.nu;
            this.dx = soln.grid.dx;
        end
        function [u_new,R_new,this] = step(this,soln,bndry_cond)
            u_new = zeros(size(soln.U,1),1);
            R_new = soln.residual(soln.U);
%             R_new = -(this.nu*(u(ip1) - 2*u(i) + u(im1))./this.dx(i).^2 ...
%                 - (u(ip1).^2 - u(im1).^2)./(4*this.dx(i)));
            u_new(this.i) = soln.U(this.i) - R_new.*soln.dt;
            u_new = bndry_cond.enforce(soln,u_new);
        end
    end
end