function plot_pos_vs_vel_rbpmf(fig,rbpmf)
figure(fig);
clf();
axis([-5 55 -35 35]);
view(-30,80);
grid on;
title('pos vs vel: x or y');
xlabel('x - position');
ylabel('vx - velocity');
hold on;

gs = Params.get_grid();

% generate gaussians to plot:
pdf_joint = zeros(gs.Nx,gs.Nvx);
x_axis = repmat(gs.x,1,gs.Nvx);
vx_axis = repmat(gs.vx,gs.Nx,1);
for x_idx = 1:size(rbpmf.xl_mean_k_k,1)
    tmp = mvnpdf(gs.vx.',rbpmf.xl_mean_k_k(x_idx),rbpmf.xl_var_k_k(x_idx));
    pdf_joint(x_idx,:) = tmp*rbpmf.Pn_k_k(x_idx);
end
plot3(x_axis.',vx_axis.',pdf_joint.');
'';
end
