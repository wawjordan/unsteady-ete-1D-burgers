classdef burgers1D_error
   properties
      i_low, i_high, imax, i, N, M
      t, t0, tf, dt {mustBeNumeric}
      order, count {mustBeNumeric}
      stencil,X,T,ptr
      x,X2,mu2
      error {mustBeNumeric}
   end
   methods
       function this = burgers1D_error( bsoln, varargin )
           
           this.i_low = bsoln.i_low;
           this.i_high = bsoln.i_high;
           this.imax = bsoln.grid.imax;
           this.count = 0;
           this.i = (bsoln.i_low:bsoln.i_high)';
           this.N = size(bsoln.grid.x,1);
           this.t0 = bsoln.t0;
           this.tf = bsoln.tf;
           this.dt = bsoln.dt;
           this.error = zeros(bsoln.grid.imax,1);
           
           defaultReconstructionOrder = 4;
           
           p = inputParser;
           validScalarPosInt = @(x) isnumeric(x) && isscalar(x) && (x > 0) && mod(x,1) == 0;
           addRequired(p,'bsoln');
           addOptional(p,'ReconstructionOrder',defaultReconstructionOrder,validScalarPosInt);
  
           parse(p,bsoln,varargin{:});
           
           this.order = p.Results.ReconstructionOrder;
           this.M = 2*ceil(this.order/2);
           this.ptr = 1:this.M+1;
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           this.x = bsoln.grid.x;
           
           this.X2 = zeros(this.M+1,this.order+1,this.imax+2);
           this.mu2 = zeros(this.imax+2,2);
           for j1 = this.i_low-1:this.i_high+1
               ind = j1-this.M/2:j1+this.M/2;
               [A,mu] = vand_matrix(this.x(ind),this.order);
               this.X2(:,:,j1) = A;
               this.mu2(j1,:) = mu;
           end
%            this.x = distributed(x);
%            this.x = x;
%            this.X2 = distributed(this.X2);
%            this.mu2 = distributed(this.mu2);
           this.X = ones(this.N,this.order+1);
           for i = 1:this.order
               this.X(:,i+1) = bsoln.grid.x.^i;
           end
           
           this.T = ones(this.M+1,this.order+1);
           this.t = zeros(this.order+1,1);
           this.stencil = zeros(this.N,this.order+1);
       end
       function [this] = reset(this)
           this.count = 0;
           this.t = zeros(this.order+1,1);
           this.stencil = zeros(this.N,this.order+1);
           this.error = zeros(this.imax,1);
           this.ptr = 1:this.M+1;
       end
       function [this] = update_time_stencil(this,bsoln)
           if bsoln.count <= this.M+1
               this.stencil(:,bsoln.count) = bsoln.U;
               this.t(bsoln.count) = bsoln.t;
           else
               this.stencil(:,this.ptr(1)) = bsoln.U;
               this.t(this.ptr(1)) = bsoln.t;
               this.ptr = circshift(this.ptr,-1);
           end
       end
%        function [u0,u1,this] = time_TE_est(this,t)
%            u0 = zeros(this.N,1);
%            u1 = zeros(this.N,1);
%            for j1 = 1:this.N
%                for j2 = 1:this.order
%                    this.T(:,j2+1) = this.t(this.ptr).^j2;
%                end
%                opts.RECT = true;
%                [P,~] = linsolve(this.T,this.stencil(j1,this.ptr)',opts);
%                P = flipud(P);
%                P1 = polyder(P);
%                u0(j1) = polyval(P,t);
%                u1(j1) = polyval(P1,t);
%            end
%        end
       function [u0,u1,this] = time_TE_est(this,t)
           u0 = zeros(this.N,1);
           u1 = zeros(this.N,1);
           [A,mu] = vand_matrix(this.t(this.ptr),this.order);
           
           for j1 = 1:this.N
               P = A\(this.stencil(j1,this.ptr)');
               pd = ddpoly(P,t,2,mu);
               u0(j1) = pd(:,1);
               u1(j1) = pd(:,2);
           end
       end
       function [u0,u1,this] = time_TE_est_mod(this,t)
%            u0 = zeros(this.N,1);
%            u1 = zeros(this.N,1);
           [U,S,V]=svd(this.stencil(:,this.ptr),0);
           [A,mu] = vand_matrix(this.t(this.ptr),this.order);
           V0 = zeros(this.M+1,1);
           V1 = zeros(this.M+1,1);
           for j1 = 1:this.M+1
               P = A\V(:,j1);
               pd = ddpoly(P,t,2,mu);
               V0(j1) = pd(:,1);
               V1(j1) = pd(:,2);
           end
           u0 = U(:,1:this.M+1)*S*V0;
           u1 = U(:,1:this.M+1)*S*V1;
       end
%        function [u0,u1,this] = time_TE_est_mod(this,t,weights)
%            u0 = zeros(this.N,1);
%            u1 = zeros(this.N,1);
%            [A,mu] = vand_matrix(this.t(this.ptr),this.order);
% %            w = eye(this.M+1);
% %            w(this.ptr(this.M),this.ptr(this.M)) = sqrt(weight);
%            weight = eye(this.M+1);
%            
%            for j1 = 1:this.N
%                weight(this.ptr(this.M+1),this.ptr(this.M+1)) = weights(j1);
%                P2 = (A'*weight*A)\(A'*weight*this.stencil(j1,this.ptr)');
% %                P1 = A\(this.stencil(j1,this.ptr)');
% %                P2 = (w*A)\(w*this.stencil(j1,this.ptr)');
% %                P2 = solve_augmented(A,this.stencil(j1,this.ptr)',weight,this.ptr(this.M+1));
%                pd = ddpoly(P2,t,2,mu);
%                u0(j1) = pd(:,1);
%                u1(j1) = pd(:,2);
%            end
%        end
       function [TE_est,R,u0] = space_TE_est(this,bsoln,u)
           u0 = zeros(this.N,1);
           u1 = zeros(this.N,1);
           u2 = zeros(this.N,1);
           for j1 = this.i_low-1:this.i_high+1
               opts.RECT = true;
               [P,~] = linsolve(this.X(j1-this.M/2:j1+this.M/2,:),...
                   u(j1-this.M/2:j1+this.M/2,:),opts);
               P = flipud(P);
               P1 = polyder(P);
               P2 = polyder(P1);
               u0(j1) = polyval(P,bsoln.grid.x(j1));
               u1(j1) = polyval(P1,bsoln.grid.x(j1));
               u2(j1) = polyval(P2,bsoln.grid.x(j1));
           end
           R = bsoln.residual(u0);
           TE_est = -bsoln.nu*u2(this.i) + u0(this.i).*u1(this.i);
       end
       function [TE_est,R,u0] = space_TE_est_mod(this,bsoln,u)
           m = this.M;
           u0 = zeros(this.N,1);
           u1 = zeros(this.N,1);
           u2 = zeros(this.N,1);
           
           for j1 = this.i_low-1:this.i_high+1
               ind = (j1-m/2:j1+m/2);
               P = this.X2(:,:,j1)\u(ind,:);
               pd = ddpoly(P,this.x(j1),3,this.mu2(j1,:));
               u0(j1) = pd(:,1);
               u1(j1) = pd(:,2);
               u2(j1) = pd(:,3);
           end
           R = bsoln.residual(u0);
           TE_est = -bsoln.nu*u2(this.i) + u0(this.i).*u1(this.i);
       end
       function [TE_est,R,u0] = space_TE_est_mod2(this,bsoln,u)
           u_mod = zeros(this.N,this.m+1);
           x_mod = zeros(this.N,this.m+1);
           for j1 = 1:this.m+1
               u_mod(:,j1) = u(this.i+j1-1-this.m/2);
               x_mod(:,j1) = this.x(this.i+j1-1-this.m/2);
           end
           u0 = zeros(this.N,1);
           u1 = zeros(this.N,1);
           u2 = zeros(this.N,1);
           
           for j1 = this.i_low-1:this.i_high+1
               ind = (j1-this.m/2:j1+this.m/2);
               P = this.X2(:,:,j1)\u(ind,:);
               pd = ddpoly(P,this.x(j1),3,this.mu2(j1,:));
               u0(j1) = pd(:,1);
               u1(j1) = pd(:,2);
               u2(j1) = pd(:,3);
           end
           R = bsoln.residual(u0);
           TE_est = -bsoln.nu*u2(this.i) + u0(this.i).*u1(this.i);
       end
%        function [TE_est,R,u0] = space_TE_est_mod2(this,bsoln,u)
%            m = this.M;
%            u0 = zeros(this.N,1);
%            u1 = zeros(this.N,1);
%            u2 = zeros(this.N,1);
%            
%            spmd
%                N1 = numlabs;
%            end
%            N1 = length(N1);
%            Numpoints = floor((this.imax+2)/N1)*ones(N1,1);
%            for j1 = 1:mod(this.imax+2,N1)
%                Numpoints(j1) = Numpoints(j1) + 1;
%            end
%            I1 = ones(N1,1);
%            I2 = ones(N1,1);
%            I1(1) = this.i_low-1;
%            I2(1) = this.i_low-1+Numpoints(1)-1;
%            for j1 = 2:N1
%                I1(j1) = I2(j1-1)+1;
%                I2(j1) = I1(j1) + Numpoints(j1)-1;
%            end
%            spmd
%                N2 = I2(labindex)-I1(labindex);
%                u02 = zeros(N2,1);
%                u12 = zeros(N2,1);
%                u22 = zeros(N2,1);
%                if labindex == 1
%                    X3 = labBroadcast(1,this.X2); % send to all
%                    mu3 = labBroadcast(1,this.mu2);
%                    x3 = labBroadcast(1,this.x);
%                    u3 = labBroadcast(1,u);
%                else
%                    X3 = labBroadcast(1);
%                    mu3 = labBroadcast(1);
%                    x3 = labBroadcast(1);
%                    u3 = labBroadcast(1);
%                end
%                
%                for j1 = I1(labindex):I2(labindex)
%                    ind = (j1-m/2:j1+m/2);
%                    ind2 = j1+1-I1(labindex);
%                    P = X3(:,:,j1)\u3(ind,:);
%                    pd = ddpoly(P,x3(j1),3,mu3(j1,:));
%                    u02(ind2) = pd(:,1);
%                    u12(ind2) = pd(:,2);
%                    u22(ind2) = pd(:,3);
%                end
%            end
%            u0(this.i_low-1:this.i_high+1) = vertcat(u02{:});
%            u1(this.i_low-1:this.i_high+1) = vertcat(u12{:});
%            u2(this.i_low-1:this.i_high+1) = vertcat(u22{:});
%            R = bsoln.residual(u0);
%            TE_est = -bsoln.nu*u2(this.i) + u0(this.i).*u1(this.i);
%        end
   end
end