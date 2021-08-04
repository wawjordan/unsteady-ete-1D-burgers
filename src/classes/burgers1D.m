classdef burgers1D < soln1D
   properties
      ExactSolutionType
      ista
      i_low, i_high, i
      error {mustBeNumeric}
      Uref, Re, nu {mustBeNumeric}
      t, t0, tf, dt {mustBeNumeric}
      CFL {mustBeNumeric}
   end
   properties (Dependent = true, SetAccess = private)
       ExactSolution
       Jacobian
       Residual
   end
   methods
       function this = burgers1D( grid, Re, varargin )
           this = this@soln1D( grid, 1 );
           this.Re = Re;
           this.Uref = 2;
           this.nu = this.Uref*this.grid.L/this.Re;
           this.i_low = this.grid.i_low;
           this.i_high = this.grid.i_high;
           this.i = (this.grid.i_low:this.grid.i_high)';
           this.error = zeros(grid.imax,1);
           
           defaultTimeAccurate = false;
           defaultCFL = 1.0;
           defaultTimeRange = [0,1];
           defaultdt = min(...
                       [abs(0.5*(this.grid.dx(this.i)./this.U(this.i))),...
                       0.499*((this.grid.dx(this.i).^2)/this.nu)],[],2);
           defaultExactSolutionType = 'steady_shock';
           expectedSolutions = {'steady_shock','unsteady_shock',...
               'pulse_plus','pulse_minus','comp_pulse'};
           
           p = inputParser;
           validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
           validScalarLogical = @(x) islogical(x) && isscalar(x);
           valid1x2 = @(x) isnumeric(x) && isrow(x) && size(x,2)==2;
           addRequired(p,'grid');
           addRequired(p,'Re',validScalarPosNum);
           addOptional(p,'TimeAccurate',defaultTimeAccurate,...
               validScalarLogical);
           addOptional(p,'CFL',defaultCFL,validScalarPosNum);
           addOptional(p,'TimeRange',defaultTimeRange,valid1x2);
           addOptional(p,'dt',defaultdt,validScalarPosNum);
           addOptional(p,'ExactSolutionType',defaultExactSolutionType,...
               @(x) any(validatestring(x,expectedSolutions)));
           
           parse(p,grid,Re,varargin{:});
           
           this.ista = p.Results.TimeAccurate;
           this.CFL = p.Results.CFL;
           this.t0 = p.Results.TimeRange(1);
           this.tf = p.Results.TimeRange(2);
           this.dt = p.Results.CFL*p.Results.dt;
           if (this.ista==true)
               this.dt = min(this.dt);
           end
           this.t = this.t0;
           
           switch(p.Results.ExactSolutionType)
               case 'steady_shock'
                   this.ExactSolutionType = @steady_shock;
               case 'unsteady_shock'
                   this.ExactSolutionType = @shock_coalesce;
               case 'pulse_plus'
                   this.ExactSolutionType = @pulse_decay_plus;
               case 'pulse_minus'
                   this.ExactSolutionType = @pulse_decay_minus;
               case 'comp_pulse'
                   this.ExactSolutionType = @compression_pulse;
           end
       end
       function ExactSolution = get.ExactSolution(this)
           ExactSolution = this.ExactSolutionType(this,this.grid.x,this.t);
       end
       function Residual = get.Residual(this)
           Residual = -(this.nu*(this.U(this.i+1) - 2*this.U(this.i) + ...
               this.U(this.i-1))./this.grid.dx(this.i).^2 - ...
               (this.U(this.i+1).^2 - this.U(this.i-1).^2)./...
               (4*this.grid.dx(this.i)));
       end
       function res = residual(this,u)
           res = -(this.nu*(u(this.i+1) - 2*u(this.i) + ...
               u(this.i-1))./this.grid.dx(this.i).^2 - ...
               (u(this.i+1).^2 - u(this.i-1).^2)./...
               (4*this.grid.dx(this.i)));
       end
       function Uex = calc_exact(this,x,t)
          Uex = this.ExactSolutionType(this,x,t);
       end
       function Jacobian = get.Jacobian(this)
           Jacobian = zeros(this.grid.imax,3);
           Jacobian(:,1) = this.nu./this.grid.dx(this.i-1).^2 + ...
               this.U(this.i-1)./(2*this.grid.dx(this.i-1));
           Jacobian(1,1) = 0;
           Jacobian(:,2) = -2*this.nu./this.grid.dx(this.i).^2;
           Jacobian(:,3) = this.nu./this.grid.dx(this.i+1).^2 - ...
               this.U(this.i+1)./(2*this.grid.dx(this.i+1));
           Jacobian(this.grid.imax,3) = 0;
       end
       function J = jacobian(this,u)
           J = zeros(this.grid.imax,3);
           J(:,1) = this.nu./this.grid.dx(this.i-1).^2 + ...
               u(this.i-1)./(2*this.grid.dx(this.i-1));
           J(1,1) = 0;
           J(:,2) = -2*this.nu./this.grid.dx(this.i).^2;
           J(:,3) = this.nu./this.grid.dx(this.i+1).^2 - ...
               u(this.i+1)./(2*this.grid.dx(this.i+1));
           J(this.grid.imax,3) = 0;
       end
       function uex = shock_coalesce(this,x,t)
           x2 = x*0.5*this.Re/this.grid.L;
           t2 = t*0.25*this.Re*this.Uref/this.grid.L;
           uex = -2*sinh(x2)./(cosh(x2)+exp(-t2));
       end
       function uex = pulse_decay_plus(this,x,t)
           x2 = x*0.5*this.Re/this.grid.L;
           t2 = t*0.25*this.Re*this.Uref/this.grid.L;
           uex = x2./t2./(1 + sqrt(t2).*exp(x2.^2./(4*t2)));
       end
       function uex = pulse_decay_minus(this,x,t)
           x2 = x*0.5*this.Re/this.grid.L;
           t2 = t*0.25*this.Re*this.Uref/this.grid.L;
           uex = x2./t2./(1 - sqrt(t2).*exp(x2.^2./(4*t2)));
       end
       function uex = steady_shock(this,x,~)
           a1 = -this.Re*this.nu/this.grid.L;
           x2 = x*0.5*this.Re/this.grid.L;
           uex = a1*tanh(x2);
       end
       function uex = compression_pulse(this,x,t)
           alpha = 2/(exp(0.5*this.Re)-1);
           z = x/(2*sqrt(t));
           uex = (2/sqrt(pi*t))*exp(-z.^2)./(alpha+erfc(z));
       end
   end
end