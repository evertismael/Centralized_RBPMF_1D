classdef PMF_bsbnd_based
   properties
       ptcls_est
       ptcls_pred
       
       X_est
       P_est
       
       X_rng
   end
   methods
      function obj = PMF_bsbnd_based()
          % get needed variables / parameters
          gs = Params.get_grid();
          [x1,x2] = ndgrid(gs.x,gs.vx);
          obj.X_rng = [x1(:),x2(:)];
          '';
      end      
      
      function obj = init(obj,X_0,P_0)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          tmp = mvnpdf(obj.X_rng,X_0.',P_0);
          obj.ptcls_est = reshape(tmp,gs.Nx,gs.Nvx);
          obj.X_est = X_0;
          obj.P_est = P_0;
          '';
          
          % how much is the norm? (normalize)
          sum_ptcls = trapz(trapz(obj.ptcls_est,1),2)*gs.dvx*gs.dx;
          obj.ptcls_est = obj.ptcls_est./sum_ptcls;
          '';
      end
      
      function obj = time_update(obj)
          % get needed variables / parameters
          gs = Params.get_grid();
          bf = Params.get_bayesian_params();
          
          % final predicted pdf:
          pred_pdf = zeros(size(obj.ptcls_est));
          for ptcl_idx = 1:size(obj.X_rng)
              mu = obj.X_rng(ptcl_idx,:);
              mu = bf.F*mu.';
              mu = mu.';
              tmp = mvnpdf(obj.X_rng,mu,bf.Q);
              p_xt_xtm1 = reshape(tmp,gs.Nx,gs.Nvx);
              
              w_tm1 = obj.ptcls_est(ptcl_idx);
              pred_pdf = pred_pdf + p_xt_xtm1*w_tm1;
              '';
          end
          obj.ptcls_pred = pred_pdf;
          '';
      end
      
      function obj = measurement_update(obj,BSs)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          % generate the likelihood grid:
          rx_diff = BSs.rx - BSs.rx_grid;
          tmp = (-0.5./BSs.vars).*(sum(conj(rx_diff).*rx_diff,5));
          log_zt_xt = sum(tmp,4);
          p_zt_xt = exp(log_zt_xt - max(log_zt_xt,[],1));
          p_zt_xt = repmat(p_zt_xt,1,gs.Nvx); % repeat along vx
          
          '';
%           figure;
%           plot(gs.x,p_zt_xt)
          
          % product particle to particle:
          obj.ptcls_est = p_zt_xt .* obj.ptcls_pred;
          
          % renormalize:
          sum_ptcls = trapz(trapz(obj.ptcls_est,1),2)*gs.dvx*gs.dx;
          obj.ptcls_est = obj.ptcls_est./sum_ptcls;
          '';
      end
      
      function obj = compute_estimates(obj)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          % means:
          x_mean = trapz(gs.x.*trapz(obj.ptcls_est,2)*gs.dvx)*gs.dx;
          vx_mean = trapz(gs.vx.*trapz(obj.ptcls_est,1)*gs.dx)*gs.dvx;
          
          obj.X_est = [x_mean;vx_mean];
          
          % vars:
          x_var = trapz((gs.x-x_mean).^2.*trapz(obj.ptcls_est,2)*gs.dvx)*gs.dx;
          vx_var = trapz((gs.vx-vx_mean).^2.*trapz(obj.ptcls_est,1)*gs.dx)*gs.dvx;
          xvx_var = trapz((gs.x-x_mean).*trapz((gs.vx-vx_mean).*obj.ptcls_est,2)*gs.dvx)*gs.dx;
          
          obj.P_est = [x_var xvx_var; xvx_var vx_var];
      end
   end
end