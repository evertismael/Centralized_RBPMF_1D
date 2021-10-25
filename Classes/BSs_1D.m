classdef BSs_1D
   properties
      tx        % tx bsbnd signal
      rx        % rx bsbnd signal
      
      rx_grid   % tx_grid (received signals as if they came from that point)
      vars      % variances at the BSs.
   end
   methods
      function obj = BSs_1D()
          %scene = Params.get_scene();
          %comm = Params.get_communication();
      end
      
      function obj = gen_tx(obj)
          bs_gs = Params.get_bs_grid();
          comm = Params.get_communication();
          
          % create Tx signal:
          bitstream = randi(2, comm.Nbps * comm.N_pilot, 1) - 1;
          obj.tx = mapping(bitstream,comm.Nbps,'qam');
          obj.tx = reshape(obj.tx,1,1,1,1,comm.N_pilot); % dim(X,Y,cfo,Nbs,Npilot)
          
          % create signal grid:
          obj.rx_grid = obj.tx.*exp(1j*comm.phi_rng.*(bs_gs.r));
          '';
      end
      
      function obj = capture_rx(obj,X_true)
          scene = Params.get_scene();
          comm = Params.get_communication();
          
          % propagation effect
          xy_target = X_true(1); % only x in this example
          delta = abs(scene.bx - xy_target); % distance from target to BSs
          delta = reshape(delta, 1,1,1,scene.Nbs);
          
          obj.rx = obj.tx.*exp(1j*comm.phi_rng.*(delta));
          
          % noise: 'more distance more noise'
          SNR_lin_ref = 10^(comm.SNR_db_ref/10);
          SNR_i = SNR_lin_ref.*((comm.delta_ref./delta).^2);
          obj.vars = reshape(1./SNR_i,1,1,1,scene.Nbs);
          w = sqrt(obj.vars/2).*(randn(size(obj.rx)) + 1j*randn(size(obj.rx)));
          
          obj.rx = obj.rx + w;
          '';
      end
      
      function X_dpe = dpe(obj)
          % get needed variables / parameters
          gs = Params.get_grid();
          
          % build the likelihood p(z|x);
          rx_diff = obj.rx - obj.rx_grid;
          tmp = (-0.5./obj.vars).*(sum(conj(rx_diff).*rx_diff,5));
          log_zt_xt = sum(tmp,4);
          p_z_xn = exp(log_zt_xt - max(log_zt_xt,[],1));
          
          p_z_xn_norm = trapz(p_z_xn,1)*gs.dx;
         
          X_dpe = trapz(gs.x.*p_z_xn,1)*gs.dx/p_z_xn_norm;
          '';
      end
   end
end