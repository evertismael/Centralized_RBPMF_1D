classdef RBPMF_pos_based
   properties
       Pn_k_km1
       Pn_k_k
       
       xl_mean_k_k    % k | k
       xl_var_k_k
       
       xl_mean_k_km1  % k | k-1
       xl_var_k_km1
       
       % Final Estimate
       X_est
       P_est
       
   end
   methods
      function obj = RBPMF_pos_based()
      end      
      
      function obj = init(obj,X_0,P_0)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          
          % X_0 = [xn;xl], P_0 = [pn pnl; pnl pl];
          % init grid in the NON-LINEAR domain:
          xn_mean_0 = X_0(1,:);
          xn_var_0 = P_0(1,1); 
          Pn_hat_0 = mvnpdf(gs.x,xn_mean_0,xn_var_0);
          c = sum(Pn_hat_0,1)*gs.dx;
          Pn_0 = Pn_hat_0./c;
          obj.Pn_k_km1 = reshape(Pn_0,gs.Nx,1);
          
                    
          % init continuous pdf per each sample in the Non-Linear space.
          xl_mean_0 = X_0(2,:);
          xl_var_0 = P_0(2,2); 
          Pn = P_0(1,1); Pnl = P_0(1,2); Pln = P_0(2,1);
          
          obj.xl_mean_k_km1 = xl_mean_0 + Pln*pinv(Pn)*(gs.x - xn_mean_0);
          obj.xl_var_k_km1 = xl_var_0 - Pln*pinv(Pn)*Pnl;
          if size(obj.xl_var_k_km1,1)==1
              obj.xl_var_k_km1 = obj.xl_var_k_km1*ones(size(obj.xl_mean_k_km1));
          end
          
          obj.X_est = X_0;
          obj.P_est = P_0;
      end
      
      function obj = time_update(obj)
          % get needed variables / parameters
          gs = Params.get_grid();
          bf = Params.get_bayesian_params();
          
          %
          Fn = bf.F(1,1); Fnl = bf.F(1,2);
          Fln = bf.F(2,1); Fl = bf.F(2,2);
          
          % VIRTUAL update
          % create virtual measurements.
          eps_j = gs.x.';
          eps_i = gs.x;
          y_ij = eps_j - Fn*eps_i; % (56)   
          c_ij = zeros(gs.Nx,gs.Nx);
          
          % normalization factor c:
          for i_idx = 1:size(gs.x,1)
            c_mean = Fnl*obj.xl_mean_k_k(i_idx);
            c_var = Fnl*obj.xl_var_k_k(i_idx)*Fnl.' + bf.Q(1,1);
            c_ij(i_idx,:) = mvnpdf(y_ij(i_idx,:).',c_mean,c_var);
          end
          
          S = Fnl*obj.xl_var_k_k*(Fnl.') + bf.Q(1,1);
          K_ij = obj.xl_var_k_k*(Fnl.')./(S);
          xl_var_k_k_st_ij = (1 - K_ij.*Fnl).*obj.xl_var_k_k;
          x1_mean_k_k_st_ij = obj.xl_mean_k_k + K_ij.*(y_ij - Fnl.*obj.xl_mean_k_k);
          '';
          % TIME update:
          xl_mean_k_km1_ij = Fl.*x1_mean_k_k_st_ij;
          xl_var_k_km1_ij = Fl.*xl_var_k_k_st_ij.*(Fl.') + bf.Q(2,2);
          '';
          
          % MOMENT MATCHING;
          P_kp1_k_ij = obj.Pn_k_k.*(c_ij).*gs.dx;
          P_hat_kp1_k_ij = P_kp1_k_ij./sum(P_kp1_k_ij,1);
          
          obj.xl_mean_k_km1 = sum(P_hat_kp1_k_ij.*xl_mean_k_km1_ij,1);
          obj.xl_var_k_km1 = sum(P_hat_kp1_k_ij.*(xl_var_k_km1_ij + (xl_mean_k_km1_ij - obj.xl_mean_k_km1).^2),1);
          obj.Pn_k_km1 = sum(P_kp1_k_ij,1);
          
          % correcting dimentions:
          obj.xl_mean_k_km1 = obj.xl_mean_k_km1.';
          obj.xl_var_k_km1 = obj.xl_var_k_km1.';
          obj.Pn_k_km1 = obj.Pn_k_km1.';
          '';
          
          % ISSUE: when P is too small the normalization factor is zero
          % SOL: look at the NaNs and reinit them.
          if sum(isnan(obj.xl_mean_k_km1),'all')>0 || sum(isnan(obj.xl_var_k_km1),'all')>0
              tmp_mean = isnan(obj.xl_mean_k_km1);
              tmp_var = isnan(obj.xl_var_k_km1);
              obj.xl_mean_k_km1(tmp_mean) = 0;
              obj.xl_var_k_km1(tmp_var) = 10;
              '';
          end
          '';
      end
      
      function obj = measurement_update(obj,Z)
          % get needed variables / parameters
          gs = Params.get_grid();
          bf = Params.get_bayesian_params();
          
          % ---NON-LINEAR part: p(z|x) grid:
          % In our case, this is a grid. In a normal case we evaluate:
          % p(z|x) = N(h(xn) + H*xl_mean, HPlH^T + R); 
          Z_rng = gs.x; % z = h(xn) + H(xl); H=0;
          p_z_xn = mvnpdf(Z_rng,Z,bf.R);
          c = sum(p_z_xn.*obj.Pn_k_km1,1)*gs.dx;
          obj.Pn_k_k = (1/c).*p_z_xn.*obj.Pn_k_km1;
          
          
          % ---LINEAR part: p(l_k_km1|n,z^k);
          % notice that in our case H=0 -> K=0;
          obj.xl_mean_k_k = obj.xl_mean_k_km1; % l_k_k = l_k_km1 + K*(z_i-z_mean);
          obj.xl_var_k_k = obj.xl_var_k_km1; %P_k_k = P_k_km1 - KHP_k_km1;
          '';
      end
      
      function obj = compute_estimates(obj)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          % LINEAR ESTIMATES: (velocities)
          xl_mean = sum(obj.Pn_k_k.*obj.xl_mean_k_k,1)*gs.dx;
          xl_var = sum(obj.Pn_k_k.*(obj.xl_var_k_k + (obj.xl_mean_k_k - xl_mean).^2),1)*gs.dx;
          
          % NONLINEAR ESTIMATES: (position)
          xn_mean = sum(obj.Pn_k_k.*gs.x,1)*gs.dx;
          xn_var = sum(obj.Pn_k_k.*(gs.x - xn_mean).^2,1)*gs.dx + gs.dx.^2/12;
          
          xnl_var = sum(obj.Pn_k_k.*(obj.xl_mean_k_k - xl_mean).*(gs.x-xn_mean),1)*gs.dx;
          
         
          obj.X_est = [xn_mean;xl_mean];
          obj.P_est = [xn_var xnl_var; xnl_var xl_var];
      end
   end
end