function target = car_1D_movement()
x0 = [10,4].';
%x0 = [10,-20].';

target = Mobile_1D(x0);
syms t

% TRAYECTORY 1: turning to right to start the roundabout.
t1 = 0; t2 = 6;%0.32;
a = 1;
w = 1;
b = 0;
law_vx = @(x_init,v_init,t_init,t) v_init -7*(cos(10*t)+ sin(0.5*10*t) - cos(0.2*10*t));
law_x = @(x_init,v_init,t_init,t) x_init + int(v_init -7*(cos(10*t)+ sin(0.5*10*t) - cos(0.2*10*t)),t_init,t);
law = Law_1D(law_x,law_vx);
target = target.add_trayectory(t1,t2,law);
end