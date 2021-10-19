classdef soln_stencil
    properties
        M, N
        max_length
        queue_length
        U, t
    end
    methods
        function this = soln_stencil( M, N, varargin )
            this.M = M;
            this.N = N;
            if nargin > 3
                error('too many input arguments')
            end
            if nargin > 2
                this.max_length = max(M,varargin{1});
            else
                this.max_length = M;
            end
            this.U = zeros(N,M);
            this.t = zeros(1,M);
            this.queue_length = 1;
        end
        function this = push(this,u_new,t_new)
            this.U(:,1:end-1) = this.U(:,2:end); % shift all columns left
            this.t(1:end-1) = this.t(:,2:end);
            this.U(:,end) = u_new; % add new vector
            this.t(end) = t_new;
            this.queue_length = min(this.M, this.queue_length+1);
        end
        function this = grow(this,k)
            if (this.M < this.max_length)
                if (this.M+k < this.max_length)
                    this.U = [zeros(this.N,k),this.U];
                    this.t = [zeros(1,k),this.t];
                    this.M = this.M+k;
                else
                    warning('specified k exceeds max_length')
                    k2 = this.max_length-this.M;
                    this.U = [zeros(this.N,k2),this.U];
                    this.t = [zeros(1,k2),this.t];
                    this.M = this.max_length;
                end
            else
                warning('stencil has max_length')
            end
        end
    end
end
            