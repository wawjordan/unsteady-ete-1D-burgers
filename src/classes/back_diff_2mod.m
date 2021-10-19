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
    
            u_new = soln.U;
            R       = nan(this.newton_max_iter,soln.neq); % allocate residual norm storage
            J       = soln.jacobian(u_new);              % calculate residual Jacobian
            res     = -soln.residual(u_new);        % calculate steady-state residual for current solution
            J2      = -(2/3)*soln.dt.*J;                  % modify jacobian for BDF2 algorithm
            J2(:,2) = J2(:,2) + 1;
            J2(1,:) = 0;
            J2(1,2) = 1; 
            J2(end,:) = 0;
            J2(end,2) = 1;
%             tmp = J2;
%             tmp(:,1) = circshift(J2(:,1),-1);
%             tmp(:,3) = circshift(J2(:,3),1);
%             J3 = spdiags(tmp,-1:1,this.imax,this.imax);
            
            % F^(k+1) = u_h^(k+1) - (4/3)*u_h^n + (1/3)*u_h^(n-1) - (2/3)dt*R_h(u_h^k)
            % RHS = -F^k
            u1    = this.um1(this.i);
            u1m1  = this.um2(this.i);
            dt = soln.dt;
            F = @(uhk) uhk - (4/3)*u1 + (1/3)*u1m1 - (2/3)*dt*res;
            R_new = -F(u_new(this.i));
            R_new(1) = soln.calc_exact(soln.grid.x(soln.i_low),soln.t)-u_new(soln.i_low);
            R_new(end) = soln.calc_exact(soln.grid.x(soln.i_high),soln.t)-u_new(soln.i_high);
            du = tridiag(J2(:,1),J2(:,2),J2(:,3),R_new);
%             du = J3\R_new;
            u_new(soln.i) = u_new(soln.i) + du;
            
            j = 1;
            if this.start == true
                this.start = false;
                this.Rnorm(1,:) = 1;
                R(1,:) = this.Rnorm;
                for k = 1:soln.neq
                    this.Rinit(1,k) = norm(R_new(:),2)/(soln.grid.imax^(1/2));
                end
            else
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(1,:) = this.Rnorm;
            end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                
                J       = soln.jacobian(u_new);              % calculate residual Jacobian
                res     = -soln.residual(u_new);              % calculate steady-state residual for current solution
                J2      = -(2/3)*soln.dt.*J;                  % modify jacobian for BDF2 algorithm
                J2(:,2) = J2(:,2) + 1;
                J2(1,:) = 0;
                J2(1,2) = 1;
                J2(end,:) = 0;
                J2(end,2) = 1;
                
                F = @(uhk) uhk - (4/3)*this.um1(this.i) + (1/3)*this.um2(this.i) - (2/3)*dt*res;
                R_new = -F(u_new(soln.i));
                R_new(1) = (soln.calc_exact(soln.grid.x(soln.i_low),soln.t)-u_new(soln.i_low));
                R_new(end) = (soln.calc_exact(soln.grid.x(soln.i_high),soln.t)-u_new(soln.i_high));
                du = tridiag(J2(:,1),J2(:,2),J2(:,3),R_new);
                
                u_new(this.i) = u_new(this.i) + du;
                
                j = j+1;
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(j,:) = this.Rnorm;
            end
            this.um2 = this.um1;

            u_new = bndry_cond.enforce(soln,u_new);
            this.um1 = u_new;
            R = R(~isnan(R(:,1)),:);
            
        end
        function res = residual(~,soln,u)
            res = 0*u(soln.i);
            ind = (soln.i_low+1:soln.i_high-1);
            res(2:soln.grid.imax-1) = -(soln.nu*(u(ind+1) - 2*u(ind) + ...
               u(ind-1))./soln.grid.dx(ind).^2 - ...
               (u(ind+1).^2 - u(ind-1).^2)./...
               (4*soln.grid.dx(ind)));
        end
        function J = jacobian(~,soln,u)
           J = zeros(soln.grid.imax,3);
           ind1 = (soln.i_low+1:soln.i_high-1);
           ind2 = (2:soln.grid.imax-1);
           J(ind2,1) = soln.nu./soln.grid.dx(ind1-1).^2 + ...
               u(ind2-1)./(2*soln.grid.dx(ind1-1));
%            J(1,1) = 0;
           J(ind2,2) = -2*soln.nu./soln.grid.dx(ind1).^2;
           J(ind2,3) = soln.nu./soln.grid.dx(ind1+1).^2 - ...
               u(ind2+1)./(2*soln.grid.dx(ind1+1));
%            J(soln.grid.imax,3) = 0;
       end
    end
end