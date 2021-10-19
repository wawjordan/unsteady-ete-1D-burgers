classdef back_diff_2_ETEmod < time_integrator_type
    properties
        imax, i_low, i_high, i
        tau
        Rnorm, Rinit
        em1, em2
        pos
    end
    methods
        function this = back_diff_2_ETEmod(soln_error)
            this.i_low = soln_error.i_low;
            this.i_high = soln_error.i_high;
            this.imax = soln_error.imax;
            this.i = soln_error.i;
            this.pos = soln_error.M;

            this.em1 = soln_error.error;
            this.em2 = soln_error.error;
%             this.u_old = soln_error.stencil;
        end
        function [e_new,R,this] = step(this,soln,soln_error,bndry_cond)
            e_new = soln_error.error;
            ind = soln_error.ptr(this.pos+1);
%             ind1 = soln_error.ptr(this.pos);
%             ind2 = soln_error.ptr(this.pos-1);
%============ use this one for asymmetric stencil     ========%
%             ind = soln_error.ptr(soln_error.M+1);
%             ind1 = soln_error.ptr(soln_error.M);
%             ind2 = soln_error.ptr(soln_error.M-1);
            
%============ use this one for symmetric stencil     ========%
%             ind = soln_error.ptr(soln_error.M/2+1);
            
%             t1 = soln_error.t(ind1);
%             t2 = soln_error.t(ind2);
%             disp(t1);
%             disp(t2);
            dt = soln.dt;
            t = soln_error.t(ind);
%             t1 = soln_error.t(ind)-dt;
%             t2 = soln_error.t(ind)-2*dt;
%             u = soln.U;
            u = soln_error.stencil(:,ind);
%             u = this.u_old(:,ind);
            
            ue = u;
            ue(this.i) = u(this.i)-e_new;
%             ue = bndry_cond.enforce(soln,ue);
            
            Rue = soln.residual(ue);
            Ru = soln.residual(u);
            
            %% LHS
            J = soln.jacobian(ue);
            J2 = -(2/3)*dt.*J;
            J2(:,2) = J2(:,2) + 1;
            J2(1,:) = 0;
            J2(1,2) = 1; 
            J2(end,:) = 0;
            J2(end,2) = 1;
            
            %% TE estimation
% % %             [u01,~,soln_error] = soln_error.time_TE_est(t1);
% % %             [u02,~,soln_error] = soln_error.time_TE_est(t2);
%             [u0,u13,soln_error] = soln_error.time_TE_est(t);
            [u0,u13,soln_error] = soln_error.time_TE_est_mod(t);
% % %             [~,~,u01] = soln_error.space_TE_est(soln,u01);
% % %             [~,~,u02] = soln_error.space_TE_est(soln,u02);
% % %             [TE,~,~] = soln_error.space_TE_est(soln,u0);
% % %             [TE,~,~] = soln_error.space_TE_est(soln,u0);
              [TE,~,~] = soln_error.space_TE_est_mod(soln,u0);
% % %             [TE,RE,u0] = soln_error.space_TE_est(soln,u0);
% % %             [~,TE] = soln_error.space_TE_est(soln,u03);
            Tdudt = u13(soln.i);
% % %             Rdudt = (3*u0(this.i) - 4*u01(this.i) + u02(this.i))/(2*dt);
% % %             Rdudt = (3*u(this.i) - 4*u1(this.i) + u2(this.i))/(2*dt);
%             u1 = soln.calc_exact(soln.grid.x,t);
%             u2 = soln.calc_exact(soln.grid.x,t1);
%             u3 = soln.calc_exact(soln.grid.x,t2);
%             TE3 = -((3*u1(this.i) - 4*u2(this.i) + u3(this.i))/(2*dt) + soln.residual(u1));
            TE = TE + Tdudt;
% % %             TE2 = RE + Rdudt;
%             clf
%             hold on;
%             plot(TE,'r')
%             plot(TE3,'k')
% % %             TE = TE-TE2;
            this.tau = TE;
            %% RHS
            RHS = -(e_new - (4/3)*this.em1 + (1/3)*this.em2 + ...
                (2/3)*dt.*(-Rue+Ru-TE));
            RHS(1) = -e_new(1);
            RHS(end) = -e_new(end);
            
            R = nan(this.newton_max_iter,soln.neq);
    
            de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            e_new = e_new + de;
            R_new = RHS;
            j = 1;
            this.Rnorm(1,:) = 1;
            R(1,:) = this.Rnorm;
            for k = 1:soln.neq
                this.Rinit(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2));
            end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                ue(this.i) = u(this.i)-e_new;
                J = soln.jacobian(ue);
                J2 = -(2/3)*dt.*J;
                J2(:,2) = J2(:,2) + 1;
                J2(1,:) = 0;
                J2(1,2) = 1;
                J2(end,:) = 0;
                J2(end,2) = 1;
                
                Rue = soln.residual(ue);
                RHS = -(e_new - (4/3)*this.em1 + (1/3)*this.em2 - ...
                    (2/3)*dt.*(Rue-Ru+TE));
                RHS(1) = -e_new(1);
                RHS(end) = -e_new(end);

                de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                e_new = e_new + de;

                R_new = RHS;
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(j+1,:) = this.Rnorm;
                j = j+1;
            end
         
            this.em2 = this.em1;
            this.em1 = e_new;
            R = R(~isnan(R(:,1)),:);
        end
        function [e_new,R,this] = step_mod(this,soln,soln_error,weight)
            e_new = soln_error.error;
            ind = soln_error.ptr(this.pos+1);
            
            dt = soln.dt;
            t = soln_error.t(ind);
            
%             u = soln.U;
            u = soln_error.stencil(:,ind);
%             u = this.u_old(:,ind);
            
            ue = u;
            ue(this.i) = u(this.i)-e_new;
%             ue = bndry_cond.enforce(soln,ue);
            
            Rue = soln.residual(ue);
            Ru = soln.residual(u);
            
            %% LHS
            J = soln.jacobian(u);
            J2 = -(2/3)*dt.*J;
            J2(:,2) = J2(:,2) + 1;
            
            %% TE estimation
            [u0,u13,soln_error] = soln_error.time_TE_est_mod(t,weight);

            [TE,~,~] = soln_error.space_TE_est(soln,u0);

            Tdudt = u13(soln.i);

            TE = TE + Tdudt;

            this.tau = TE;
            %% RHS
            RHS = -(e_new - (4/3)*this.em1 + (1/3)*this.em2 + ...
                    (2/3)*dt.*(-Rue+Ru-TE));
            
            R = nan(this.newton_max_iter,soln.neq);
            
            de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            e_new = e_new + de;
            R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
            j = 1;
                this.Rnorm(1,:) = 1;
                R(1,:) = this.Rnorm;
                for k = 1:soln.neq
                    this.Rinit(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2));
                end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                ue(this.i) = u(this.i)-e_new;

                J = soln.jacobian(ue);
                J2 = -(2/3)*dt.*J;
                J2(:,2) = J2(:,2) + 1;
                
                Rue = soln.residual(ue);
                RHS = -(e_new - (4/3)*this.em1 + (1/3)*this.em2 - ...
                    (2/3)*dt.*(Rue-Ru+TE));
                de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                e_new = e_new + de;
                R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(j,:) = this.Rnorm;
                j = j+1;
            end
            this.em2 = this.em1;
            this.em1 = e_new;
            R = R(~isnan(R(:,1)),:);
        end
    end
end