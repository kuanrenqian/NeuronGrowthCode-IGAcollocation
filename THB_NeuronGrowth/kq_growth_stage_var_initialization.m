% iter_stage2_begin = 1000;
iter_stage2_begin = 2000;
iter_stage3_begin = 8000;
iter_stage45_begin = 15000;
rot_iter_invl = 500;
% phi_actin = reshape(cm.NuNv*phi,lenu,lenv);
% param = GetParam(phi_actin,dtime);
% actin_start = rot_iter_invl*2;

Rot = zeros(1,20);
rotate = zeros(1,20);
rotate_intv = zeros(1,20);
size_Max = 0;

Max_x = 0;
Max_y = 0;
delta_L = 1;
term_change = 1;

var_save_invl = 500;
png_save_invl = 100;