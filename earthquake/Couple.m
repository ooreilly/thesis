classdef Couple

  properties

    % Objects for each side
    S1; 
    S2;

    % Penalty parameters
    a;
    b;

    % Fault and friction parameters
    fault;
      
    % Enable debugging, 
    % plots the two functions involved in the nonlinear solve
    use_debug = false;

  end

  methods

    % s1 is object to the left and s2 is the object to the right
    function obj = Couple(S1,S2,fault)

      obj.S1 = S1;
      obj.S2 = S2;
      obj.fault = fault;

      % Penalty parameters
      obj.a = obj.S1.rho*obj.S1.c;
      obj.b = obj.S2.rho*obj.S2.c;

    end

    % Couples two sides assuming the interface is locked
    function [vs ss] = couple_locked(obj,u1,u2)
      
      % Get v
      v1 = obj.S1.get_v(u1);
      v1 = v1(end);
      v2 = obj.S2.get_v(u2);
      v2 = v2(1);

      % Get sigma
      s1 = obj.S1.get_s(u1);
      s1 = s1(end);
      s2 = obj.S2.get_s(u2);
      s2 = s2(1);

      vs =           1/(obj.a+obj.b)*(s2 - s1) + 1/(obj.a+obj.b)*(obj.a*v2 + obj.b*v1);
      ss = obj.a*obj.b/(obj.a+obj.b)*(v2 - v1) + 1/(obj.a+obj.b)*(obj.b*s2 + obj.a*s1);

    end

    % Couples two sides with fault dynamics
    function [vs1 vs2 ss] = couple(obj,u1,u2,Psi,t)

      % Get v
      v1 = obj.S1.get_v(u1);
      v1 = v1(end);
      v2 = obj.S2.get_v(u2);
      v2 = v2(1);

      % Get sigma
      s1 = obj.S1.get_s(u1);
      s1 = s1(end);
      s2 = obj.S2.get_s(u2);
      s2 = s2(1);

      a = obj.a; 
      b = obj.b;

      f = obj.fault;

      V   = v2 - v1;
      c   = a*b/(a+b);
      phi = (a*s2+b*s1)/(a+b);
      d   = phi + c*V;

      % Sigma^*
      ss = @(Vs) d - c*Vs;

      g  = @(Vs) ss(Vs) - f.F(Vs,Psi) + f.taul(t) + f.tau0;

      bracket  =  d/c+(f.taul(t)+f.tau0)/c;
      lb = min([0,bracket]);
      ub = max([0,bracket]);

      obj.debug(c,d,@f.F,ss,Psi);
      Vs = secant_method(g,V,lb,ub,1e-6,1e-9,obj.use_debug);

      % Set shear stresses
      ss  = f.F(Vs,Psi)-f.tau0;

      % Set velocities
      vs1 = v1 - (s1 - ss)/a;
      vs2 = v2 + (s2 - ss)/b;

    end

    function V = slip_velocity(obj,u1,u2)
      % Get v
      v1 = obj.S1.get_v(u1);
      v1 = v1(end);
      v2 = obj.S2.get_v(u2);
      v2 = v2(1);

      V = v2 - v1;

    end

    function debug(obj,c,d,F,ss,Psi)
        if ~obj.use_debug
          return;
        end

        v = linspace(-abs(d/c),abs(d/c),100);
        plot(v,F(v,Psi),'o-',v,ss(v));
        drawnow;
        legend('F(V)','aV+b');
    end

  end










end
