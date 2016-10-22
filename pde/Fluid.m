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
      obj.sizes = [n+1,n+2];

      obj = obj.interior();
      obj = obj.mechanical_energy();
      obj = obj.sat_self();

    end

    function obj = material_properties(rho,G)

      obj.rho = rho;
      obj.G   = G;
      obj.c   = sqrt(rho/G);
      obj.Z   = obj.rho*obj.c;

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
      Ln     = -obj.mu/obj.rho*obj.ppi*r.e0p*r.e0m'*obj.Dm*n;
      Ld     = -obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.e0m*r.e0p'*n;

      % Right side
      n     =  1;
      Rn     = -obj.mu/obj.rho*obj.ppi*r.eNp*r.eNm'*obj.Dm*n;
      Rd     = -obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.eNm*r.eNp'*n;


      obj.A = obj.A + Ln + Ld + Rn + Rd;

    end

    function obj = neumann_left(obj)
      
      r     = restrictions(obj.n);
      n     = -1;
      Ld    = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.e0m*r.e0p'*n;
      obj.A = obj.A + Ld;
      

    end

    function obj = neumann_right(obj)
      
      r     = restrictions(obj.n);
      n     = 1;
      Rd    = obj.mu/obj.rho*inv(obj.Pp)*obj.Dm'*r.eNm*r.eNp'*n;
      obj.A = obj.A + Rd;
      

    end
    
    function obj = absorbing_bc_right(obj)
      
      r = restrictions(obj.n);
      n = 1;
      
      % Add sigma^*
      S     = -0.5*obj.Z*obj.ppi/obj.rho*r.eNp*r.eNp';
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,1,1,S);
      S     =  0.5*obj.pmi/obj.rho*r.eNp*r.eNm'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,1,2,S);

      % Add v^*
      S     = 0.5*obj.ppi/obj.rho*r.eNm*r.eNp'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,2,1,S);
      S     = -0.5*obj.pmi/obj.rho/obj.Z*r.eNm*r.eNm';
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,2,2,S);


    end

    % The vector returned by either coupling function is added to the right hand side 
    % of the system to enforce any valid coupling conditions
    % For example, 
    % du1/dt = A*u1 + couple_right(vs1,ss1)
    % du2/dt = B*u2 + couple_left(vs2,ss2)
    function u = couple_left(obj,vs,ss)
      n                  = -1;
      u                  = zeros(sum(obj.sizes),1);
      u(1)               =   obj.ppi/obj.rho*ss*n;
      u(obj.sizes(1)+1)  =   obj.pmi*obj.G*vs*n;
    end
    
    function u = couple_right(obj,vs,ss)
      n                  = 1;
      u                  = zeros(sum(obj.sizes),1);
      u(obj.sizes(1))    =   obj.ppi/obj.rho*ss*n;
      u(end)             =   obj.pmi*obj.G*vs*n;
    end

    function obj = mechanical_energy(obj)

      [xp,xm,obj.Pp,obj.Pm,Qp,Qm] = sbp_staggered(obj.order,obj.n,obj.h);
      obj.H = 0.5*obj.rho*obj.Pp;

    end

    function info = test_stability(obj,halt)
      [info.is_stable, info.is_energy_stable,info.eig_s,info.eig_es] = ...
      test_energy_stability(obj.A,obj.H,1,halt,1e-8);
    end

    function v_ = get_v(obj,u)
      v_ = u(1:obj.n+1);
    end
    
    function s_ = get_s(obj,u)
      s_ = u(obj.n+2:end);
    end





  end







end
