classdef soln_type
   properties
      L {mustBeNumeric}
      N {mustBeNumeric}
      dx {mustBeNumeric}
      x {mustBeNumeric}
      u {mustBeNumeric}
      Re {mustBeNumeric}
      Uref {mustBeNumeric}
      nu {mustBeNumeric}
      t {mustBeNumeric}
      t0 {mustBeNumeric}
      tf {mustBeNumeric}
      dt {mustBeNumeric}
      BC {mustBeNumeric}
      CFL {mustBeNumeric}
      R {mustBeNumeric}
      J {mustBeNumeric}
      time_integrator
      residuals
      eps = 1e-6;
      maxiter = 5000;
      maxnewtoniter = 10;
   end
   methods
       function obj = Soln_Type(L,N,Re,Uref)
           obj.L = L;
           obj.N = N;
           obj.Re = Re;
           obj.Uref = Uref;
           obj.nu = Uref*L/Re;
           obj.dx = L/(N-1);
           obj.x = (-L/2:L/(N-1):L/2)';
       end
       function obj = setup_t_unsteady(obj,t0,tf,dt)
           obj.t0 = t0;
           obj.tf = tf;
           obj.dt = dt;
           obj.t = (t0:dt:tf)';
       end
       function obj = setup_t_steady(obj)
           obj.t0 = 0;
           obj.tf = 0;
           obj.dt = 0;
           obj.t = zeros(obj.maxiter,1);
       end
       function uex = shock_coalesce(obj,x,t)
           x2 = x*0.5*obj.Re/obj.L;
           t2 = t*0.25*obj.Re*obj.Uref/obj.L;
           uex = -2*sinh(x2)./(cosh(x2)+exp(-t2));
       end
       function uex = pulse_decay_plus(obj,x,t)
           x2 = x*0.5*obj.Re/obj.L;
           t2 = t*0.25*obj.Re*obj.Uref/obj.L;
           uex = x2./t2./(1 + sqrt(t2).*exp(x2.^2./(4*t2)));
       end
       function uex = pulse_decay_minus(obj,x,t)
           x2 = x*0.5*obj.Re/obj.L;
           t2 = t*0.25*obj.Re*obj.Uref/obj.L;
           uex = x2./t2./(1 - sqrt(t2).*exp(x2.^2./(4*t2)));
       end
       function uex = steady_shock(obj,x)
           a1 = -obj.Re*obj.nu/obj.L;
           x2 = x*0.5*obj.Re/obj.L;
           uex = a1*tanh(x2);
       end
   end
end