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
pmf = PMF_bsbnd_based();
%pmf = PMF_bsbnd_based_v2();

% Base stations:
BSs = BSs_1D();
BSs = BSs.gen_tx();


% params and history
N_t = size(target.t_vect,2);

X_true_hist = zeros(2,N_t);
X_est_hist = zeros(2,N_t);
P_est_hist = zeros(4,N_t);

X_dpe_hist = zeros(1,N_t); %  just position

% figures
f_pos = figure('Position', [4 73 560 922]);
f_vel_pos = figure('Position', [568 72 560 922]);
f_vel_pos_pdf = figure('Position', [1130 72 560 922]);

comm = Params.get_communication();
init = Params.get_initials();


for t_idx = 1:N_t
    X_true = target.history(:,t_idx)
    
    % create measurement:
    BSs = BSs.capture_rx(X_true);
        
    if t_idx == 1
        X_0 = X_true;
        P_0 = init.P_0;
        pmf = pmf.init(X_0,P_0);
        pmf.ptcls_pred = pmf.ptcls_est;
    else
        % correct
        pmf = pmf.measurement_update(BSs);  
        pmf = pmf.compute_estimates();
        
        % predict
        pmf = pmf.time_update();
    end
   
    
    % save history:
    X_true_hist(:,t_idx) = X_true;
    X_dpe_hist(:,t_idx) = BSs.dpe();
    
    X_est_hist(:,t_idx) = pmf.X_est;  
    P_est_hist(:,t_idx) = pmf.P_est(:);
    
    % plot:
    if mod(t_idx,5)==0
        plot_state_vector(f_pos,target.t_vect(:,1:t_idx),...
            X_true_hist(:,1:t_idx), X_est_hist(:,1:t_idx),...
            X_dpe_hist(:,1:t_idx), P_est_hist(:,1:t_idx));
        
        plot_pos_vs_vel(f_vel_pos,X_true_hist(:,1:t_idx),...
            X_est_hist(:,1:t_idx));  
        plot_pdf_pos_vs_vel(f_vel_pos_pdf,pmf);  
        pause(0.1);
    end

end





