function plot_rbpmf_info(fig,rbpmf)
figure(fig);
clf();
subplot(2,2,1);

gs = Params.get_grid();
% Pn_k_k
plot(gs.x,rbpmf.Pn_k_k,'DisplayName','k|k'); hold on;
plot(gs.x,rbpmf.Pn_k_km1,'DisplayName','k|k-1');
xlim([-5 55]);
title('Pn');
xlabel('x');
legend();

% means xl
subplot(2,2,2);
plot(gs.x,rbpmf.xl_mean_k_k,'DisplayName','k|k');hold on;
plot(gs.x,rbpmf.xl_mean_k_km1,'DisplayName','k|k-1');
xlim([-5 55]);
ylim([-300, 300]);
title('xl-means');
%legend();

% variances xl
subplot(2,2,3);
plot(gs.x,rbpmf.xl_var_k_k,'DisplayName','k|k'); hold on;
plot(gs.x,rbpmf.xl_var_k_km1,'DisplayName','k|k-1');
xlim([-5 55]);
ylim([0, 100]);
title('xl-vars');
legend();
'';
end
