function plot_pos_vs_vel(fig,X_true_hist,X_est_hist)
figure(fig);
clf();
axis([-5 55 -35 35]);
grid on;
title('position vs velocity');
hold on;

plot(X_true_hist(1,:),X_true_hist(2,:),'DisplayName','true');
plot(X_est_hist(1,:),X_est_hist(2,:),'DisplayName','hat');
legend();
end