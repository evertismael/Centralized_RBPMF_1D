classdef PMF_bsbnd_based_v2
   properties
       ptcls_est
       ptcls_pred
       
       X_est
       P_est
       
       X_rng
   end
   methods
      function obj = PMF_bsbnd_based_v2()
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
          bf = Params.get_bayesian_params();
          gs = Params.get_grid();

          % final predicted pdf:
          X_t_tm1 = bf.F*obj.X_est;
          P_t_tm1 = bf.F*obj.P_est*bf.F.' + bf.Q;

          tmp = mvnpdf(obj.X_rng,X_t_tm1.',P_t_tm1);
          obj.ptcls_pred = reshape(tmp,gs.Nx,gs.Nvx);
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