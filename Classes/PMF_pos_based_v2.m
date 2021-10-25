classdef PMF_pos_based_v2
    % This class does the same of the original PMF v1, except that it does
    % the time update in an analytical way.
   properties
       ptcls_est
       ptcls_pred
       
       X_est
       P_est

       X_rng
   end
   methods
      function obj = PMF_pos_based_v2()
          % get needed variables / parameters
          gs = Params.get_grid();
          [x1,x2] = ndgrid(gs.x,gs.vx);
          obj.X_rng = [x1(:),x2(:)];
          '';
      end      
      
      function obj = init(obj,X_0,P_0)
          % get needed variables / parameters
          gs = Params.get_grid();
          bf = Params.get_bayesian_params();
          
          tmp = mvnpdf(obj.X_rng,X_0.',P_0);
          obj.ptcls_est = reshape(tmp,gs.Nx,gs.Nvx);
          obj.X_est = X_0;
          obj.P_est = P_0;
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
      
      function obj = measurement_update(obj,Y)
          % get needed variables / parameters
          gs = Params.get_grid();
          bf = Params.get_bayesian_params();
          
          Y_rng = obj.X_rng(:,1); % y = H*x
          tmp = mvnpdf(Y_rng,Y,bf.R);
          
          % generate the likelihood grid:
          p_yt_xt =  reshape(tmp,gs.Nx,gs.Nvx);
          
          % product particle to particle:
          obj.ptcls_est = p_yt_xt .* obj.ptcls_pred;
          
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