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

    % Fluid 
    f;
    n1;

    % Solid
    s;
    n2;

    % Coupling
    C;

    % Save pgf data
    % name: figure name (name_field_num)
    % path: figure path
    %  num: counter for current snapshot in time
    name = 'fluid-solid'; 
    path = 'pgf_data';
    num  = 1;


  end

  methods

  
    function obj = Solver(f,s)

      obj = obj.initialize(f,s);
      
    end

    function obj = solve(obj,CFL,tend,show_plot,frame_stride,save_at_time)

      % Time step and compute number of time steps to take
      obj.dt = CFL*min([obj.f.cfl() obj.s.cfl()]);
      obj.nt = round(tend/obj.dt);

      obj.frame_stride = frame_stride;
      f = @(x,t) obj.problem(x,t);
      t = 0;
      for i = 1:obj.nt
        obj.u = lsrk4(f,obj.u,t,obj.dt);
        %obj   = save_fields(obj,t);
        %obj   = save_data(obj,t,save_at_time);
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
        u2 = obj.u(obj.n1+1:end);
        plot(obj.f.xp-1,u1,obj.s.xp,obj.s.get_v(u2));
        drawnow;

    end

    function obj = initialize(obj,f,s)

      % Fluid
      f     = f.initialize();
      obj.f = f;
      obj.f = obj.f.neumann_left();
      %obj.f = obj.f.neumann_right();

      % Solid
      s     = s.initialize();
      obj.s = s;
      
      obj.n1 = sum(obj.f.sizes);
      obj.n2 = sum(obj.s.sizes);
    
      % Boundary conditions
      %obj.s = obj.s.absorbing_bc_left();
      obj.s = obj.s.absorbing_bc_right();

      % Initialize coupling
      obj.C = Couple(obj.f,obj.s);

      % Initial conditions
      xp = obj.f.xp;
      a = 0.1;
      u1 = exp(-(xp-0.5).^2/(2*a^2))';

      xp = obj.s.xp; xm = obj.s.xm;
      u2 = [exp(-(xp-0.5).^2/(2*a^2))'; 0*xm'];
      obj.u = [u1;u2];

      obj.n1 = length(u1);
      obj.n2 = length(u2);

    end


    % Function that sets the right hand side for each time step
    function u = problem(obj,u,t)

      u1 = u(1:obj.n1); 
      u2 = u(obj.n1+1:end);
      u  = [obj.f.A*u1;obj.s.A*u2];
      %% Couple
      [vs ss]      = obj.C.couple(u1,u2);
      u1c          = obj.f.couple_right(vs,ss,u1,t);
      u2c          = obj.s.couple_left(vs,ss);
      uc = [u1c;u2c];
      u  = u + uc;


    end

    function obj = save_fields(obj,t)
      
      if ~obj.save_now(t);
        return;
      end

      u1  = obj.u(1:obj.n1); 
      u2  = obj.u(obj.n1+1:end-1);

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
