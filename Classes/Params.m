classdef Params
    methods (Static)
        
        function glb = get_globals()
            glb.T = (1e-1);     % sampling Time
            glb.c  = 2.998e8;   % speed of light
        end
        
        function init = get_initials()
            init = {};
            init.P_0 = 15*eye(2);
            init.P_0(2,2) = 100;
        end
                
        function bf = get_bayesian_params()
            glb = Params.get_globals();
            
            bf = {};
            % noise variance
            var_v = 20^2; % max acceleration change in the interval
            bf.Q = [0.5*(glb.T^2); glb.T]*var_v*([0.5*(glb.T^2); glb.T].');
            tmp = .1*randn(2,2); % to make it invertible;
            bf.Q = bf.Q + tmp*tmp.';
            
            % process model:
            glb = Params.get_globals();
            bf.F = [1 glb.T; 0 1];
            
            % measurement model:
            bf.R = 2*eye(1);
        end
            
        function comm = get_communication()
            glb = Params.get_globals();
            
            comm = {};
            comm.fc = 2e9;
            comm.B = 40e6;
            comm.N_subcrr = 1024;
            comm.pilot_spacing = 16;
            comm.N_pilot = comm.N_subcrr/comm.pilot_spacing;
            comm.deltaPhi = 2*pi*comm.B/(comm.N_pilot*glb.c);
            comm.Nbps = 2;
            comm.pilot_idx_rng = (1:comm.N_pilot);
            comm.phi_rng = comm.deltaPhi*comm.pilot_idx_rng;
            comm.phi_rng = reshape(comm.phi_rng,[1,1,1,1,comm.N_pilot]);
            
            
            comm.SNR_db_ref = -10; % reference SNR
            comm.delta_ref = 20; % reference distance for SNR
            
            % noise variances for z = h(x) + wn;
            comm.xy_var_n = (.1)^2; % when h(x) = H*x; H = [1 0]
            % comm.var_n = 1; % when h(x) = received signal% We use a
            % reference SNR in dB
            
        end
        
        function gs = get_grid()
          gs = {};
          gs.x_min = 0; gs.x_max = 50; gs.Nx = 100;
          gs.vx_min = -30; gs.vx_max = 30; gs.Nvx = 50;
          
          gs.dx = (gs.x_max - gs.x_min)/gs.Nx;
          gs.dvx = (gs.vx_max - gs.vx_min)/gs.Nvx;
          
          gs.x = gs.x_min:gs.dx:gs.x_max-gs.dx;
          gs.vx = gs.vx_min:gs.dvx:gs.vx_max-gs.dvx;
          
          % reshape grid:
          gs.x = reshape(gs.x,gs.Nx,1);
          gs.vx = reshape(gs.vx,1,gs.Nvx);
        end
        
        function bs_gs = get_bs_grid()
            scene = Params.get_scene();
            gs = Params.get_grid();
            
            % distance grid to each BSs
            bs_gs.r = abs(gs.x - reshape(scene.bx,1,1,1,scene.Nbs)); 
        end
            
        function scene = get_scene()
            scene.bx = [0 0 50 50];  % positions of BSs
            scene.Nbs = size(scene.bx,2);
        end
    end
end