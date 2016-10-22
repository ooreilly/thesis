classdef Fluid

  properties
    n; 
    h;
    sizes;
    A;
    H;
    order;

    % Grids
    xp; xm;

    % Quadrature rules
    Pp; Pm;
    % Difference operators
    Dp; Dm;
    % Corner weights
    ppi, pmi;

    % Material properties
    rho = 1, mu = 1;

  end

  methods

    function obj = Fluid(n,order)

      obj.order = order;
      obj.n     = n;
      obj.h     = 1/n;
      obj.sizes = n+1;

    end

    function obj = initialize(obj)

      obj = obj.interior();
      obj = obj.mechanical_energy();
      obj = obj.sat_self();

    end

    function obj = material_properties(rho,mu)

      obj.rho = rho;
      obj.mu  = mu;

    end

    function obj = interior(obj)

      [obj.xp,obj.xm,Pp,Pm,Qp,Qm] = sbp_staggered(obj.order,obj.n,obj.h);
      obj.Dp = inv(Pp)*Qp; obj.Dm = inv(Pm)*Qm;
      obj.Pp = Pp; obj.Pm = Pm;

      obj.ppi = 1/obj.Pp(1,1);
      obj.pmi = 1/obj.Pm(1,1);

      obj.A = obj.mu/obj.rho*obj.Dp*obj.Dm;

    end

    % Estimates the timestep using the CFL condition
    function dt = cfl(obj)

      dt = obj.h^2/(obj.mu/obj.rho);

    end

    function obj = sat_self(obj)

      r = restrictions(obj.n);

      % Left side
      n     = -1;
      Ln    = -obj.mu/obj.rho*obj.ppi*r.e0p*r.e0m'*obj.Dm*n;
      Ld    = -obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.e0m*r.e0p'*n;

      % Right side
      n     =  1;
      Rn    = -obj.mu/obj.rho*obj.ppi*r.eNp*r.eNm'*obj.Dm*n;
      Rd    = -obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.eNm*r.eNp'*n;

      obj.A = obj.A + Ln + Ld + Rn + Rd;

    end

    function obj = neumann_left(obj)
      
      r     = restrictions(obj.n);
      n     = -1;
      L    = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.e0m*r.e0p'*n;
      obj.A = obj.A + L;

    end

    function obj = neumann_right(obj)
      
      r     = restrictions(obj.n);
      n     = 1;
      R    = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.eNm*r.eNp'*n;
      obj.A = obj.A + R;

    end
    
    function obj = dirichlet_left(obj)
      
      r     = restrictions(obj.n);
      n     = -1;
      L     = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.e0m*r.e0p'*n;
      obj.A = obj.A + L;

    end

    function obj = dirichlet_right(obj)
      
      r     = restrictions(obj.n);
      n     = 1;
      R     = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.eNm*r.eNp'*n;
      obj.A = obj.A + R;

    end

    % The vector returned by either coupling function is added to the right hand side 
    % of the system to enforce any valid coupling conditions
    % For example, 
    % du1/dt = A*u1 + couple_right(us1,taus1)
    % du2/dt = B*u2 + couple_left(us2,taus2)
    function u = couple_left(obj,us,taus)

      r     = restrictions(obj.n);
      n     = -1;
      ctau  = 1/obj.rho*obj.ppi*r.e0p*n*taus;
      cu    = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.e0m*n*us;
      u     = ctau + cu;

    end
    
    function u = couple_right(obj,us,taus,uu,t)

      r     = restrictions(obj.n);
      n     = 1;
      ctau  = 1/obj.rho*obj.ppi*r.eNp*n*taus;
      cu    = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.eNm*n*us;
      u     = ctau + cu;

    end

    function s = sat_parameter(obj)

      s = obj.mu/obj.h;

    end

    function obj = mechanical_energy(obj)

      [xp,xm,obj.Pp,obj.Pm,Qp,Qm] = sbp_staggered(obj.order,obj.n,obj.h);
      obj.H = 0.5*obj.rho*obj.Pp;

    end

    function info = test_stability(obj,halt)

      [info.is_stable, info.is_energy_stable,info.eig_s,info.eig_es] = ...
      test_energy_stability(obj.A,obj.H,1,halt,1e-8);

    end

    function tau = get_tau(obj,u)

      tau = obj.mu*obj.Dm*u;

    end

  end

end
