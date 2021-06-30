classdef trapezoid_method_ETE < time_integrator_type
    properties
        imax, i_low, i_high, i
        tau
        Rnorm, Rinit
        pos, u_old
    end
    methods
        function this = trapezoid_method_ETE(soln_error)
            this.i_low = soln_error.i_low;
            this.i_high = soln_error.i_high;
            this.imax = soln_error.imax;
            this.i = soln_error.i;
            this.pos = soln_error.M;
            this.u_old = soln_error.stencil;
        end
        function this = set.u_old(this,value)
            this.u_old = value;
        end
        function [e_new,R,this] = step(this,soln,soln_error,~)
%============ use this one for asymmetric stencil     ========%
%             ind1 = soln_error.ptr(soln_error.M);
%             ind2 = soln_error.ptr(soln_error.M+1);            
%============ use this one for symmetric time stencil ========%
%             ind1 = soln_error.ptr(soln_error.M/2);
%             ind2 = soln_error.ptr(soln_error.M/2+1);

            ind1 = soln_error.ptr(this.pos);
            ind2 = soln_error.ptr(this.pos+1);
            
            t1 = soln_error.t(ind1);
            t2 = soln_error.t(ind2);
            
            u1 = soln_error.stencil(:,ind1);
            u2 = soln_error.stencil(:,ind2);
            
            ue1 = u1;
            ue2 = u2;
            ue1(this.i) = u1(this.i)-soln_error.error;
            ue2(this.i) = u2(this.i)-soln_error.error;
            
            Rue1 = soln.residual(ue1);
            Rue2 = soln.residual(ue2);
            
            Ru1 = soln.residual(u1);
            Ru2 = soln.residual(u2);
            
            [u01,u11] = soln_error.time_TE_est(t1);
            [u02,u12] = soln_error.time_TE_est(t2);
            
            [TE1,~] = soln_error.space_TE_est(soln,u01);
            [TE2,~] = soln_error.space_TE_est(soln,u02);
            
            TE1 = TE1 + u11(this.i);
            TE2 = TE2 + u12(this.i);
            
            F1 = Rue1 - Ru1 + TE1;
            F2 = Rue2 - Ru2 + TE2;
            
%             [TE1,RE1,u01] = soln_error.space_TE_est(soln,u01);
%             [TE2,RE2,u02] = soln_error.space_TE_est(soln,u02);
%             Tdudt1 = u11(this.i);
%             Tdudt2 = u12(this.i);
%             Rdudt = (u02(this.i)-u01(this.i))/soln.dt;
%             TE1 = TE1 + Tdudt1;
%             TE2 = TE2 + Tdudt2;
%             TE1 = TE1 + (RE1+Rdudt);
%             TE2 = TE2 + (RE2+Rdudt);
%             F1 = Rue1 - Ru1 + TE1;
%             F2 = Rue2 - Ru2 + TE2;
            
            this.tau = 0.5*(TE1+TE2);
            
            J = soln.jacobian(u1);
            J2 = -(1/2)*soln.dt.*J;
            J2(:,2) = J2(:,2) + 1;
            
            R = nan(this.newton_max_iter,soln.neq);

            RHS = 0.5*soln.dt*(F1+F2);
            
            de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            
            e_new = soln_error.error + de;
            
            R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
            
            j = 1;
            this.Rnorm(1,:) = 1;
            R(1,:) = this.Rnorm;
            for k = 1:soln.neq
                this.Rinit(1,k) = norm(R_new(:,k));
            end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                
                J = soln.jacobian(ue2);
                J2 = -(1/2)*soln.dt.*J;
                J2(:,2) = J2(:,2) + 1;
                
                ue2(this.i) = u2(this.i)-e_new;
                
                Rue2 = soln.residual(ue2);
                
                F2 = Rue2 - Ru2 + TE2;
                
                RHS = 0.5*soln.dt*(F1+F2)-(e_new-soln_error.error);
                
                de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                
                e_new = e_new + de;
                
                R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
                
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/this.Rinit(1,k);
                end
                R(j,:) = this.Rnorm;
                j = j+1;
            end
            
%             TE3 = 0.5*(soln.residual(soln.calc_exact(soln.grid.x,t1)) + soln.residual(soln.calc_exact(soln.grid.x,t2)) );
%             TE3 = TE3 + (soln.calc_exact(soln.grid.x(soln.i),t2)-soln.calc_exact(soln.grid.x(soln.i),t1))/soln.dt;
%             clf;
%             hold on;
%             plot(soln.grid.x(soln.i),-TE3,'k');
%             plot(soln.grid.x(soln.i),0.5*(TE1+TE2),'g');
%             hold off;
            R = R(~isnan(R(:,1)),:);
        end
    end
end