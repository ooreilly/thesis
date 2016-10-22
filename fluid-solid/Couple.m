classdef Couple

  properties

    % Objects for each side
    S1; 
    S2;

    % Penalty parameters
    a;
    b;
      
    % Enable debugging, 
    use_debug = false;

  end

  methods

    % s1 is object to the left and s2 is the object to the right
    function obj = Couple(S1,S2)

      obj.S1 = S1;
      obj.S2 = S2;

      % Penalty parameters
      obj.a = obj.S1.sat_parameter();
      obj.b = obj.S2.sat_parameter();

    end

    % Couples the fluid on the left to the solid on the right
    function [vs ss] = couple(obj,fluid,solid)
      
      % Get v
      u   = fluid(end);
      tau = obj.S1.get_tau(fluid);
      tau = tau(end);

      % Get sigma
      v = obj.S2.get_v(solid);
      v = v(1);
      s = obj.S2.get_s(solid);
      s = s(1);

      vs =           1/(obj.a+obj.b)*(s - tau) + 1/(obj.a+obj.b)*(obj.b*v  + obj.a*u);
      ss = obj.a*obj.b/(obj.a+obj.b)*(v - u)   + 1/(obj.a+obj.b)*(obj.a*s  + obj.b*tau);

    end

  end










end
