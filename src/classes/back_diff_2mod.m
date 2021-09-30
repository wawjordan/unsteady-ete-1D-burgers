classdef back_diff_2mod < time_integrator_type
    properties
        imax, i_low, i_high, i
        Rnorm, Rinit
        um1, um2
        start
    end
    methods
        function this = back_diff_2mod(soln)
            this.i_low = soln.i_low;
            this.i_high = soln.i_high;
            this.imax = soln.grid.imax;
            this.i = soln.i;
            this.um1 = soln.calc_exact(soln.grid.x,soln.t - min(soln.dt));
            this.um2 = soln.calc_exact(soln.grid.x,soln.t - 2*min(soln.dt));
            this.start = true;
        end
        function [u_new,R,this] = step(this,soln,bndry_cond)
    
            u_new = zeros(size(soln.U,1),soln.neq);
            R = nan(this.newton_max_iter,soln.neq);
            J = soln.jacobian(soln.U);
            res = soln.residual(soln.U);
            J2 = -(2/3)*soln.dt.*J;
            J2(:,2) = J2(:,2) + 1;
            RHS = -soln.U(this.i) + (4/3)*this.um1(this.i) - ...
                (1/3)*this.um2(this.i) - (2/3)*soln.dt.*res;
            
            %%% Attempt at incorporating extrapolation BC
%             a1 = -(2/3)*soln.dt*(...
%                 soln.nu/(soln.grid.dx(this.i_low)^2) + ...
%                 (2*soln.U(this.i_low) - soln.U(this.i_low+1))/...
%                 (soln.grid.dx(this.i_low) + soln.grid.dx(this.i_low+1)));
%             cN = -(2/3)*soln.dt*(...
%                 soln.nu/(soln.grid.dx(this.i_high)^2) + ...
%                 (2*soln.U(this.i_high) - soln.U(this.i_high-1))/...
%                 (soln.grid.dx(this.i_high) + soln.grid.dx(this.i_high-1)));
%             % Left boundary
%             J2(1,2) = J2(1,2) + 2*a1;
%             J2(1,3) = J2(1,3) - a1;
%             % Right boundary
%             J2(end,1) = J2(end,1) - cN;
%             J2(end,2) = J2(end,2) + 2*cN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            soln.U = bndry_cond.enforce(soln,soln.U);
            a1 = -(2/3)*soln.dt*(...
                soln.nu/(soln.grid.dx(this.i_low)^2) + ...
                soln.U(this.i_low-1)/...
                (soln.grid.dx(this.i_low) + soln.grid.dx(this.i_low+1)));
            cN = -(2/3)*soln.dt*(...
                soln.nu/(soln.grid.dx(this.i_high)^2) + ...
                soln.U(this.i_high+1)/...
                (soln.grid.dx(this.i_high) + soln.grid.dx(this.i_high-1)));
            RHS(1) = RHS(1) - a1*(soln.U(this.i_low-1)  - this.um1(this.i_low-1));
            RHS(end) = RHS(end) - cN*(soln.U(this.i_high+1) - this.um1(this.i_high+1));
            
            
            
            du = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            u_new(this.i) = soln.U(this.i) + du;
            R_new = RHS - (J2(:,1).*du + J2(:,2).*du + J2(:,3).*du);
            %%
            u_new = bndry_cond.enforce(soln,u_new);
            j = 1;
            if this.start == true
                this.start = false;
                this.Rnorm(1,:) = 1;
                R(1,:) = this.Rnorm;
                for k = 1:soln.neq
                    this.Rinit(1,k) = norm(R_new(:,k))/(soln.grid.imax^(1/2));
                end
            else
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(1,:) = this.Rnorm;
            end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                res = soln.residual(u_new);
                J = soln.jacobian(u_new);
                J2 = -(2/3)*soln.dt.*J;
                J2(:,2) = J2(:,2) + 1;
                RHS = -u_new(this.i) + (4/3)*this.um1(this.i) - ...
                (1/3)*this.um2(this.i) - (2/3)*soln.dt.*res;
            
                u_new = bndry_cond.enforce(soln,u_new);
                a1 = -(2/3)*soln.dt*(...
                    soln.nu/(soln.grid.dx(this.i_low)^2) + ...
                    u_new(this.i_low-1)/...
                    (soln.grid.dx(this.i_low) + soln.grid.dx(this.i_low+1)));
                cN = -(2/3)*soln.dt*(...
                    soln.nu/(soln.grid.dx(this.i_high)^2) + ...
                    u_new(this.i_high+1)/...
                    (soln.grid.dx(this.i_high) + soln.grid.dx(this.i_high-1)));
                RHS(1) = RHS(1) - a1*(soln.U(this.i_low-1)  - this.um1(this.i_low-1));
                RHS(end) = RHS(end) - cN*(soln.U(this.i_high+1) - this.um1(this.i_high+1));

                du = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                u_new(this.i) = u_new(this.i) + du;
                R_new = RHS - (J2(:,1).*du + J2(:,2).*du + J2(:,3).*du);
                %%
%                 u_new = bndry_cond.enforce(soln,u_new);
                j = j+1;
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(j,:) = this.Rnorm;
            end
            this.um2 = this.um1;
            this.um1 = u_new;
            R = R(~isnan(R(:,1)),:);
        end
    end
end