classdef grid1D < handle
    properties
        x, dx {mustBeNumeric}
        n_ghost, i_low, i_high {mustBeNumeric}
        xmin, xmax {mustBeNumeric}
        imax {mustBeNumeric}
        L {mustBeNumeric}
    end
    methods
        function this = grid1D( x, n_ghost )
            this.xmin = min(x);
            this.xmax = max(x);
            this.imax = length(x);
            this.L = max(x)-min(x);
            if nargin < 2
                this.n_ghost = 0;
                this.i_low = 1;
                this.i_high = this.imax;
                this.x = x(:);
                xplus = [(2*x(1)-x(2));x(:);(2*x(end)-x(end-1))];
                this.dx = 0.5*(diff(xplus(1:end-1))+diff(xplus(2:end)));
            else
                this.n_ghost = n_ghost;
                this.i_low = 1+n_ghost;
                this.i_high = this.imax+n_ghost;
                dxlow = x(1)-x(2);
                dxhigh = x(end)-x(end-1);
                xlow = flipud((x(1)+dxlow:dxlow:x(1)+(n_ghost)*dxlow)');
                xhigh = (x(end)+dxhigh:dxhigh:x(end)+(n_ghost)*dxhigh)';
                this.x = [xlow;x(:);xhigh];
                xplus = [(2*this.x(1)-this.x(2));this.x(:);...
                    (2*this.x(end)-this.x(end-1))];
                this.dx = 0.5*(diff(xplus(1:end-1))+diff(xplus(2:end)));
            end
        end
    end
end