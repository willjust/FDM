D_1 = 8/12; D_2 = 1; D_3 = 15/12; 
L_1 = 1500; L_2 = 800; L_3 = 1200;
v_2 = 20*4/pi / (sqrt((f_2/f_1) * (L_2 / L_1) * (D_1/ D_2))*D_1^2 + D_2^2 + sqrt((f_3/f_2) * (L_3/L_2) * (D_3/D_2)) * D_3^2);
v_1 = sqrt(f_2/f_1 * L_2/L_1 * D_1 / D_2 * v_2^2);
v_3 = sqrt(f_2/f_3 * L_2/L_3 * D_3 / D_2 * v_2^2);
fprintf('V1 = %3.3f \t V2= %3.3f \t V3 = %3.3f\n', v_1, v_2, v_3);
Re_1 = v_1 * D_1 / 2.0919E-5 * 1.94;
Re_2 = v_2 * D_2 / 2.0919E-5 * 1.94;
Re_3 = v_3 * D_3 / 2.0919E-5 * 1.94;
f_1 = moody(0.04/8, Re_1); f_2 = moody(0.04/12, Re_2);  f_3 = moody(0.04/15, Re_2);
DP1 = f_1 * L_1 /D_1 * 1.94 / 2 * v_1^2;
DP2 = f_2 * L_2 /D_2 * 1.94 / 2 * v_2^2;
DP3 = f_3 * L_3 /D_3 * 1.94 / 2 * v_3^2;

%% Q3. 
Q = 0.06684; D = 0.824/12; k = 5.88E-5; e = k/D; L_1 = 20; L_2 = 40; g = 32.2;
v_1 = Q * 4 / pi / D^2 / (1 + sqrt(f_1/f_2 * L_1/L_2 + 538.5996 *pi / 16 * g * D^3 / (f_2*1.94*L_2)));
v_2 = sqrt(f_1/f_2 * L_1/L_2 + 538.5996*g*D^3*pi/16 / (f_2*1.94*L_2))*v_1;
Re_1 = v_1 * D / 2.0919E-5 * 1.94; 
Re_2 = v_2 * D / 2.0919E-5 * 1.94;
f_1 = moody(e, Re_1);
f_2 = moody(e, Re_2);
DP1 = f_1*L_1/D * 1.94/ 2 * v_1^2;
DP2 = f_2*L_2/D * 1.94/ 2 * v_2^2;
Q_1 = v_1 * pi * D^2/4 / 0.002228; Q_2 = v_2*pi*D^2/4 / 0.002228;
fprintf('V1 = %3.4f \t V2= %3.4f\n', v_1, v_2);
fprintf('DP1 = %3.4f \t DP2= %3.4f\n', DP1, DP2);
fprintf('Q1 = %3.4f \t Q2= %3.4f\n', Q_1, Q_2);