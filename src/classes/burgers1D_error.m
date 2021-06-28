classdef burgers1D_error
   properties
      i_low, i_high, imax, i, N, M
      t, t0, tf, dt {mustBeNumeric}
      order, count {mustBeNumeric}
      stencil,X,T,ptr
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
       function [u0,u1,this] = time_TE_est(this,t)
           u0 = zeros(this.N,1);
           u1 = zeros(this.N,1);
           for j1 = 1:this.N
               for j2 = 1:this.order
                   this.T(:,j2+1) = this.t(this.ptr).^j2;
               end
               opts.RECT = true;
               [P,~] = linsolve(this.T,this.stencil(j1,this.ptr)',opts);
               P = flipud(P);
               P1 = polyder(P);
               u0(j1) = polyval(P,t);
               u1(j1) = polyval(P1,t);
           end
       end
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
   end
end