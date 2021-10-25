classdef RBPMF_bsbnd_based_v2
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
      function obj = RBPMF_bsbnd_based_v2()
          '';
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
          bf = Params.get_bayesian_params();
                    
          % compute mean and variances:
          obj = obj.compute_estimates();
          
          % predict based on the moments:
          X_k_km1 = bf.F*obj.X_est;
          P_k_km1 = bf.F*obj.P_est*(bf.F.') + bf.Q;
          
          obj = obj.init(X_k_km1,P_k_km1);
          '';
      end
      
      function obj = measurement_update(obj,BSs)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          % ---NON-LINEAR part: p(z|x) grid:
          % In our case, this is a grid. In a normal case we evaluate:
          % p(z|x) = N(h(xn) + H*xl_mean, HPlH^T + R); 
          
          % *** generate the likelihood grid (for my case):
          rx_diff = BSs.rx - BSs.rx_grid;
          tmp = (-0.5./BSs.vars).*(sum(conj(rx_diff).*rx_diff,5));
          log_zt_xt = sum(tmp,4);
          p_z_xn = exp(log_zt_xt - max(log_zt_xt,[],1));
          
          % proceed as before:
          c = trapz(p_z_xn.*obj.Pn_k_km1,1)*gs.dx;
          obj.Pn_k_k = (1/c).*p_z_xn.*obj.Pn_k_km1;
          
          
          % ---LINEAR part: p(l_k_km1|n,z^k);
          % notice that in our case H=0 -> K=0;
          obj.xl_mean_k_k = obj.xl_mean_k_km1; % l_k_k = l_k_km1 + K*(z_i-z_mean);
          obj.xl_var_k_k = obj.xl_var_k_km1; %P_k_k = P_k_km1 - KHP_k_km1;
          
          % --- JOINT distribution:
          % product of gaussians xL and the sample weithgs Pn_k_k;
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