classdef Fault

      properties

        % Rate-and-state friction
        %   f0: Static friction coefficient
        %   V0: Reference sliding velocity
        %    a: Direct effect parameter
        %    b: Evolution effect parameter
        %    L: State-evolution distance
        % Psi0: Initial state
        f0   = 0.6; 
        V0   = 1e-6;
        a    = 0.01;
        b    = 0.015;
        L    = 0.1;
        Psi0 = 0.52;

        % Fault properties
        % tau0: Background shear stress
        % taul: Fault load function
        %       e.g., taul = @(t) 1-exp(-t);
        tau0  = 0.2;
        taul;

      end

      methods

        function obj = Fault()

        end

        % Steady state friction
        function y = fss(obj,V)
          if abs(V)< eps
            y = obj.f0;
            return;
          end
          y = obj.f0 - (obj.b - obj.a)*sign(V)*log(abs(V)/obj.V0);
        end

        % Friction law
        function y = F(obj,V,Psi)
          y = obj.a*asinh(V/(2*obj.V0)*exp(Psi/obj.a));
        end

        % Evolution law
        function y = G(obj,V,Psi)
          y = -V/obj.L*(obj.F(V,Psi) - obj.fss(V));
        end


      end



end
