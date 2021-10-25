function plot_state_vector(fig,t_vect,X_true_hist,X_est_hist,Z_hist, P_est_hist)
figure(fig);
clf();

subplot(3,1,1)
title('position: x vs t');
ylim([-5 55]);
xlim([0 6]);
grid on;
hold on;
plot(t_vect,(X_true_hist(1,:)),'DisplayName','true');
plot(t_vect,X_est_hist(1,:),'DisplayName','hat')
plot(t_vect,Z_hist,'.','DisplayName','Z/DPE');
legend();

subplot(3,1,2)
title('velocity: vx vs t');
ylim([-30 30]);
xlim([0 6]);
grid on;
hold on;
plot(t_vect,X_true_hist(2,:),'DisplayName','true');
plot(t_vect,X_est_hist(2,:),'DisplayName','hat');
legend();


subplot(3,1,3)
title('STD: sqrt(P) vs t');
ylim([0 15]);
xlim([0 6]);
grid on;
hold on;
plot(t_vect,sqrt(P_est_hist));
legend({'x','xVx','Vxx','vx'});
end