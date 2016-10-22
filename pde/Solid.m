classdef Solid


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
    % Corner weights
    ppi, pmi;

    % Material properties
    rho = 1, G = 1; c = 1; Z = 1;

  end

  methods

    function obj = Solid(n,order)

      obj.order = order;
      obj.n     = n;
      obj.h     = 1/n;
      obj.sizes = [n+1,n+2];

    end

    function obj = initialize(obj)

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

    % Estimates the time step using the CFL condition
    function dt = cfl(obj)

      dt = obj.h/obj.c;

    end

    function obj = interior(obj)

      [obj.xp,obj.xm,Pp,Pm,Qp,Qm] = sbp_staggered(obj.order,obj.n,obj.h);
      Dp = inv(Pp)*Qp; Dm = inv(Pm)*Qm;

      obj.A = block_matrix(obj.sizes,obj.sizes);
      obj.A = block_matrix_insert(obj.A,obj.sizes,obj.sizes,1,2,1/obj.rho*Dp);
      obj.A = block_matrix_insert(obj.A,obj.sizes,obj.sizes,2,1,obj.G*Dm);

    end

    function obj = sat_self(obj)

      r = restrictions(obj.n);

      obj.ppi = 1/obj.Pp(1,1);
      obj.pmi = 1/obj.Pm(1,1);

      % Left side
      n     = -1;
      S     = -obj.ppi/obj.rho*r.e0p*r.e0m'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,1,2,S);
      S     = -obj.pmi*obj.G*r.e0m*r.e0p'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,2,1,S);

      % Right side
      n     = 1;
      S     = -obj.ppi/obj.rho*r.eNp*r.eNm'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,1,2,S);
      S     = -obj.pmi*obj.G*r.eNm*r.eNp'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,2,1,S);

    end
    
    function s = sat_parameter(obj)

      s = obj.Z;

    end

    function obj = absorbing_bc_left(obj)
      
      r = restrictions(obj.n);
      n = -1;

      % Add sigma^*
      S     = -0.5*obj.Z*obj.ppi/obj.rho*r.e0p*r.e0p';
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,1,1,S);
      S     =  0.5*obj.pmi/obj.rho*r.e0p*r.e0m'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,1,2,S);

      % Add v^*
      S     = 0.5*obj.ppi/obj.rho*r.e0m*r.e0p'*n;
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,2,1,S);
      S     = -0.5*obj.pmi/obj.rho/obj.Z*r.e0m*r.e0m';
      obj.A = block_matrix_add(obj.A,obj.sizes,obj.sizes,2,2,S);

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
      obj.H = block_matrix(obj.sizes,obj.sizes);
      obj.H = block_matrix_insert(obj.H,obj.sizes,obj.sizes,1,1,0.5*obj.rho*obj.Pp);
      obj.H = block_matrix_insert(obj.H,obj.sizes,obj.sizes,2,2,0.5/obj.G*obj.Pm);

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
