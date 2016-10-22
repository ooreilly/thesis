classdef Solver

  properties

    % Time stepping options and properties
    % frame_stride: how often to plot or save (integer)
    CFL;
    dt;
    nt;
    frame_stride;

    % Solution
    u;

    % Solid 1 
    s1;
    n1;

    % Solid 2
    s2;
    n2;

    % Coupling
    C;
    friction;

    % Fault variables as a function of time
    % tau: fault shear strength
    % Psi: rate and state friction variable
    tau;
    Psi;
    t;

    % Save pgf data
    % name: figure name (name_field_num)
    % path: figure path
    %  num: counter for current snapshot in time
    name = 'earthquake'; 
    path = 'pgf_data';
    num  = 1;


  end

  methods

  
    function obj = Solver(CFL,n1,n2,order,friction,tend)

      obj.friction = friction;
      obj = obj.initialize(n1,n2,order);
      
      % Time step and compute number of time steps to take
      obj.dt = CFL*min([obj.s1.h obj.s2.h]);
      obj.nt = round(tend/obj.dt);


    end

    function obj = solve(obj,show_plot,frame_stride,save_at_time)

      obj.frame_stride = frame_stride;
      f = @(x,t) obj.problem(x,t);
      t = 0;
      for i = 1:obj.nt
        obj.u = lsrk4(f,obj.u,t,obj.dt);
        obj   = save_fields(obj,t);
        obj   = save_data(obj,t,save_at_time);
          t   = i*obj.dt;

        plot_solution(obj,show_plot,t);

        info(obj,i,t,frame_stride);

      end

    end

    function info(obj,i,t,frame_stride)
        if mod(i,frame_stride) ~= 0
          return;
        end

        num_pad = 1 + num2str(round(log10(obj.nt)));
        fprintf(['i = %0' num_pad 'd \t t = %.2f \t %.2f %% \n'],...
                i,t,double(i)/double(obj.nt)*100);

    end

    function plot_solution(obj,show_plot,t)
        if ~show_plot
          return
        end
        
        if ~obj.save_now(t);
          return;
        end

        u1 = obj.u(1:obj.n1); 
        u2 = obj.u(obj.n1+1:end-1);
        plot(obj.s1.xp-1,obj.s1.get_v(u1),obj.s1.xp,obj.s1.get_v(u2));
        drawnow;
    end

    function obj = initialize(obj,n1,n2,order)

      obj.n1 = n1;
      obj.n2 = n2;

      % Discretize each side of the fault
      obj.s1 = Solid(obj.n1,order);
      obj.s1 = obj.s1.initialize();
      obj.s2 = Solid(obj.n2,order);
      obj.s2 = obj.s2.initialize();
    
      % Boundary conditions
      obj.s1 = obj.s1.absorbing_bc_left();
      obj.s2 = obj.s2.absorbing_bc_right();

      % Initialize coupling
      obj.C = Couple(obj.s1,obj.s2,obj.friction);

      % Initial conditions
      xp = obj.s1.xp; xm = obj.s1.xm;
      a = 0.1;
      u1 = [0*exp(-(xp-0.5).^2/(2*a^2))'; 0*xm'];
      xp = obj.s2.xp; xm = obj.s2.xm;
      u2 = [0*exp(-(xp-0.5).^2/(2*a^2))'; 0*xm'];
      obj.u = [u1;u2;obj.friction.Psi0];
      obj.n1 = length(u1);
      obj.n2 = length(u2);

      % Initialize storage fault variables
      tau = [];
      Psi = [];

    end


    % Function that sets the right hand side for each time step
    function u = problem(obj,u,t)

      u1 = u(1:obj.n1); 
      u2 = u(obj.n1+1:end-1);
      Psi = u(end);
      V  = obj.C.slip_velocity(u1,u2);
      u  = [obj.s1.A*u1;obj.s2.A*u2;obj.friction.G(V,Psi)];
      % Couple
      [vs1 vs2 ss] = obj.C.couple(u1,u2,Psi,t);
      u1c          = obj.s1.couple_right(vs1,ss);
      u2c          = obj.s2.couple_left(vs2,ss);
      uc = [u1c;u2c;Psi];
      u  = u + uc;


    end

    function obj = save_fields(obj,t)
      
      if ~obj.save_now(t);
        return;
      end

      u1  = obj.u(1:obj.n1); 
      u2  = obj.u(obj.n1+1:end-1);
      Psi = obj.u(end);
      [vs1 vs2 ss] = obj.C.couple(u1,u2,Psi,t);
      obj.tau(end+1) = ss;
      obj.Psi(end+1) = Psi;
      obj.t          = t;

    end
    
    function obj = save_data(obj,t,save_at_times)
      if ~any(abs(t-save_at_times) < obj.dt/2)
        return
      end

      u1  = obj.u(1:obj.n1); 
      u2  = obj.u(obj.n1+1:end-1);
      Psi = obj.u(end);

      pgf.x1  = obj.s1.xp-1;
      pgf.v1  = obj.s1.get_v(u1);
      pgf.x2  = obj.s2.xp;
      pgf.v2  = obj.s2.get_v(u2);

      filename = @(field) sprintf('%s_%s_%d',obj.name,field,obj.num);

      export_pgf(pgf,obj.path,filename('v'));

      obj.num = obj.num + 1;

    end
    
    function b = save_now(obj,t)
      b = false;
      if mod(round(t/obj.dt),obj.frame_stride) == 0;
        b = true;
      end
    end

  end

end
