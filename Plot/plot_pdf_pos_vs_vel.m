function plot_pdf_pos_vs_vel(fig,pmf)
figure(fig);
clf();
axis([-5 55 -35 35]);
grid on;
title('pos vs vel: x or y');
hold on;
tmp = pmf.ptcls_est;
imagesc([0,49],[-30,29],fliplr(rot90(squeeze(tmp),-1)));
end