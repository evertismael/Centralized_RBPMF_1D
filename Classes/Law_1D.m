classdef Law_1D
   properties
      law_x
      law_vx
   end
   methods
      function obj = Law_1D(law_x,law_vx)
          obj.law_x = law_x;
          obj.law_vx = law_vx;
      end
   end
end