clear all; clc; close all;
addpath('Classes')  
addpath('Targets')
addpath('Plot')  
addpath('Libs') 
rng(1);
% -------------------------------------------------------------------------
% 1.- Define & Generate target trajetories
% -------------------------------------------------------------------------
glb = Params.get_globals();
if 1==1
    target = car_up_left_v3();
    target = target.gen_trayectory(glb.T);
    % adapt to 1D:
    target.x0 = target.x0([1,2],:);
    target.history = target.history([1,2],:);
else
    target = car_1D_movement();
    target = target.gen_trayectory(glb.T);
end

% trackers:
rbpmf = RBPMF_pos_based();
%rbpmf = RBPMF_pos_based_v2();

% params and history
N_t = size(target.t_vect,2);

X_true_hist = zeros(2,N_t);
X_est_hist = zeros(2,N_t);
P_est_hist = zeros(4,N_t);
Z_hist = zeros(1,N_t);
% figures
f_pos = figure('Position', [4 73 560 922]);
f_vel_pos = figure('Position', [568 72 560 922]);
f_vel_pos_pdf = figure('Position', [1130 573 560 420]);
f_info = figure('Position', [1132 65 560 420]);

comm = Params.get_communication();
init = Params.get_initials();

for t_idx = 1:N_t
    X_true = target.history(:,t_idx)
    Z = X_true(1);
    Z = Z + sqrt(comm.xy_var_n)*randn(size(Z));
    
    if t_idx == 1
        X_0 = X_true;
        P_0 = init.P_0;
        rbpmf = rbpmf.init(X_0,P_0);
        
        % Both est and pred grids are initiated with the same values.
        rbpmf.Pn_k_k = rbpmf.Pn_k_km1;
        rbpmf.xl_mean_k_k = rbpmf.xl_mean_k_km1;
        rbpmf.xl_var_k_k = rbpmf.xl_var_k_km1;
        '';
    else
        % correct
        rbpmf = rbpmf.measurement_update(Z);
        % predict
        rbpmf = rbpmf.time_update();
        rbpmf = rbpmf.compute_estimates();
    end
   
    
    % save history:
    X_true_hist(:,t_idx) = X_true;
    Z_hist(:,t_idx) = Z;
    
    X_est_hist(:,t_idx) = rbpmf.X_est;  
    P_est_hist(:,t_idx) = rbpmf.P_est(:);
    
    % plot:
    if mod(t_idx,1)==0
        plot_state_vector(f_pos,target.t_vect(:,1:t_idx),...
            X_true_hist(:,1:t_idx), X_est_hist(:,1:t_idx),...
            Z_hist(:,1:t_idx), P_est_hist(:,1:t_idx));
        plot_pos_vs_vel(f_vel_pos,X_true_hist(:,1:t_idx),...
            X_est_hist(:,1:t_idx));  
        
        plot_pos_vs_vel_rbpmf(f_vel_pos_pdf,rbpmf)
        
        plot_rbpmf_info(f_info,rbpmf);
        pause(0.1);
    end

end
