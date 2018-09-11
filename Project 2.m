%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Multivariable Control
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Develop several (at least 2) reduced order models
% (a) Explore several model reduction techniques
% (b) Compare wrt. differing inputs (step, impulse, etc.)
% (c) Determine which models work best with given inputs
% (d) Determine which models work best for a broad range of frequencies (Bode or frequency
% response analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

A = [0      1     0      0      0     0      0      0     0      0      0      0;
    -0.202 -1.150 0      0      0     0      0      0     0      0      0      0;
     0      0     0      1      0     0      0      0     0      0      0      0;
     0      0     0      0      1     0      0      0     0      0      0      0;
     0      0    -2.360 -13.60 -12.80 0      0      0     0      0      0      0;
     0      0     0      0      0     0      1      0     0      0      0      0;
     0      0     0      0      0     0      0      1     0      0      0      0;
     0      0     0      0      0    -1.620 -9.400 -9.150 0      0      0      0;
     0      0     0      0      0     0      0      0     0      1      0      0;
     0      0     0      0      0     0      0      0     0      0      1      0;
     0      0     0      0      0     0      0      0     0      0      0      1;
     0      0     0      0      0     0      0      0    -188.0 -111.6 -116.4 -20.80];
 
B = [0      0;
     1.0439 4.1486;
     0      0; 
     0      0; 
    -1.794  2.6775;
     0      0;
     0      0;
     1.0439 4.1486;
     0      0;
     0      0; 
     0      0;
    -1.704 2.6775];
 
C = [0.2640 0.8060 -1.420 -15.00 0 0      0      0      0      0      0      0;
     0      0       0      0     0 4.9000 2.1200 1.9500 9.3500 25.800 7.1400 0];

D = [0 0;
     0 0]; 

[Model,Diagonal] = eig(A); 
 
% poles identification
sys = ss(A,B,C,D)
[p,z] = pzmap(sys);

%% ------- Residue Approach Model Reduction ------- %%

% Make Diagonal Matrix in Descending Order
M_hat = [Model(:,1),Model(:,4),Model(:,7),Model(:,9),Model(:,10),Model(:,8),Model(:,5),Model(:,2),Model(:,6),Model(:,12),Model(:,11),Model(:,3)];
 
Ad = inv(M_hat)*A*M_hat;
Bd = inv(M_hat)*B;
Cd = C*M_hat;

%%------- Model 1 (Fifth Order System) -------%%
%{
A11_1 = Ad(1:5,1:5);
A12_1 = Ad(1:5,6:12);
A21_1 = Ad(6:12,1:5);
A22_1 = Ad(6:12,6:12);

B1_1  = Bd(1:5,:);
B2_1  = Bd(6:12,:);

C1_1  = Cd(:,1:5);
C2_1  = Cd(:,6:12);

% Calculating Reduced Order Matrices
Ar_1  = A11_1 - A12_1*inv(A22_1)*A21_1;
Br_1  = B1_1 - A12_1*inv(A22_1)*B2_1;
Cr_1  = C1_1 - C2_1*inv(A22_1)*A21_1;
Dr_1  = D - C2_1*inv(A22_1)*B2_1;

% State-Space Representation of Reduced Order System
sys_r1 = ss(Ar_1,Br_1,Cr_1,Dr_1)

% Using direct function reduce()
dir1 = reduce(sys,5);

figure;
step(sys);
hold on
step(sys_r1);
title('Step Response of Original and Reduced 5th Order');

figure;
impulse(sys);
hold on
impulse(sys_r1);
title('Impulse Response of Original and Reduced 5th Order');

figure;
step(sys-sys_r1);
title('Step Response Difference between Original and Reduced 5th Order');

figure;
impulse(sys-sys_r1);
title('Impulse Response Difference between Original and Reduced 5th Order');

figure;
step(dir1);
hold on
step(sys_r1);
title('Step Response using Direct Function and Calculated Reduced 5th Order');

% Frequency Response 
figure;
bode(sys,sys_r1);
title('Frequency Response of Original and Reduced System 1');

%%------- Model 2 (Sixth Order System) -------%%

A11_2 = Ad(1:6,1:6);
A12_2 = Ad(1:6,7:12);
A21_2 = Ad(7:12,1:6);
A22_2 = Ad(7:12,7:12);

B1_2  = Bd(1:6,:);
B2_2  = Bd(7:12,:);

C1_2  = Cd(:,1:6);
C2_2  = Cd(:,7:12);

% Calculating Reduced order matrices
Ar_2  = A11_2 - A12_2*inv(A22_2)*A21_2;
Br_2  = B1_2 - A12_2*inv(A22_2)*B2_2;
Cr_2  = C1_2 - C2_2*inv(A22_2)*A21_2;
Dr_2  = D - C2_2*inv(A22_2)*B2_2;

% State-Space Representation of Reduced Order System
sys_r2 = ss(Ar_2,Br_2,Cr_2,Dr_2)

% Using direct function reduce()
dir2 = reduce(sys,6);

figure;
step(sys);
hold on
step(sys_r2);
title('Step Response of Original and Reduced 6th Order');

figure;
impulse(sys);
hold on
impulse(sys_r2);
title('Impulse Response of Original and Reduced 6th Order');

figure;
step(sys-sys_r2);
title('Step Response Difference between Original and Reduced 6th Order');

figure;
impulse(sys-sys_r2);
title('Impulse Response Difference between Original and Reduced 6th Order');

figure;
step(dir2);
hold on
step(sys_r2);
title('Step Response using Direct Function and Calculated Reduced 6th Order');

% Frequency Response 
figure;
bode(sys,sys_r2);
title('Frequency Response of Original and Reduced System 2');

%%------- Model 3 (Seventh Order System) -------%%

A11_3 = Ad(1:7,1:7);
A12_3 = Ad(1:7,8:12);
A21_3 = Ad(8:12,1:7);
A22_3 = Ad(8:12,8:12);

B1_3 = Bd(1:7,:);
B2_3 = Bd(8:12,:);

C1_3 = Cd(:,1:7);
C2_3 = Cd(:,8:12);

% Calculating Reduced order matrices
Ar_3 = A11_3 - A12_3*inv(A22_3)*A21_3;
Br_3 = B1_3 - A12_3*inv(A22_3)*B2_3;
Cr_3 = C1_3 - C2_3*inv(A22_3)*A21_3;
Dr_3 = D - C2_3*inv(A22_3)*B2_3;

% State-Space Representation of Reduced Order System
sys_r3 = ss(Ar_3,Br_3,Cr_3,Dr_3)

% Using direct function reduce()
dir3 = reduce(sys,7);

figure;
step(sys);
hold on
step(sys_r3);
title('Step Response of Original and Reduced 7th Order');

figure;
impulse(sys);
hold on
impulse(sys_r3);
title('Impulse Response of Original and Reduced 7th Order');

figure;
step(sys-sys_r3);
title('Step Response Difference between Original and Reduced 7th Order');

figure;
impulse(sys-sys_r3);
title('Impulse Response Difference between Original and Reduced 7th Order');

figure;
step(dir3);
hold on
step(sys_r3);
title('Step Response using Direct Function and Calculated Reduced 7th Order');

% Frequency Response 
figure;
bode(sys,sys_r3);
title('Frequency Response of Original and Reduced System 3');
%}
%% ------- Continued Fraction Model Reduction ------- %%
%%------- Second Order System -------%%

% Controller Type Block Companion Form
[n,m] = size(B);
disp(n);
disp(m);
gamma = n/m;        % should be integer value

% Controllability Matrix
P = [B,A*B,A^2*B,A^3*B,A^4*B,A^5*B];
inv_P = inv(P);     % exists

Bc = [zeros(2);zeros(2);zeros(2);zeros(2);zeros(2);eye(2)];
Tc1 = Bc'*inv_P;
Tc = [Tc1; Tc1*A; Tc1*A^2; Tc1*A^3; Tc1*A^4; Tc1*A^5];

Ac = Tc*A*inv_P;
Bc = Tc*B;
Cc = C*inv_P;
   
% State-Space Representation for Controller Type Block Companion Form
sys_bc = ss(Ac,Bc,Cc,0);

% First row of Routh Matrix Algorithm is from Dr(s)
A11 = -Ac(11:12,1:2);
A12 = -Ac(11:12,3:4);
A13 = -Ac(11:12,5:6);
A14 = -Ac(11:12,7:8);
A15 = -Ac(11:12,9:10);
A16 = -Ac(11:12,11:12);
A17 = eye(2);

% Second row of Routh Matrix Algorithm is from Nr(s)
A21 = Cc(1:2,1:2);
A22 = Cc(1:2,3:4);
A23 = Cc(1:2,5:6);
A24 = Cc(1:2,7:8);
A25 = Cc(1:2,9:10);
A26 = Cc(1:2,11:12);
 
% Calculations of Hankel Coefficients
H1  = A11*inv(A21)      % H1
A31 = A12-(H1*A22);
A32 = A13-(H1*A23);
A33 = A14-(H1*A24);
A34 = A15-(H1*A25);
A35 = A16-(H1*A26);

H2  = A21*inv(A31)      % H2
A41 = A22-(H2*A32);
A42 = A23-(H2*A33);
A43 = A24-(H2*A34);
A44 = A25-(H2*A35);

H3  = A31*inv(A41)      % H3
A51 = A32*(H3*A42);
A52 = A33*(H3*A43);
A53 = A34*(H3*A44);

H4  = A41*inv(A51)      % H4
A61 = A42*inv(A52);
A62 = A43*inv(A53);

% Calculating Reduced order matrices 
A_con = [-H1*H2 -H1*H4; -H1*H2 -(H1+H3)*H4];
B_con = [eye(2);eye(2)];
C_con = [H2 H4];

% State-Space Representation of Reduced Order System
sys_con = ss(A_con,B_con,C_con,0)

figure;
step(sys);
hold on
step(sys_con);
title('Step Response of Original and Continued Fraction Reduced Model');

figure;
impulse(sys);
hold on
impulse(sys_con);
title('Impulse Response of Original and Continued Fraction Reduced Model');

figure;
step(sys-sys_con);
title('Step Response Difference between Original and Continued Fraction Reduced Model');

figure;
impulse(sys-sys_con);
title('Impulse Response Difference between Original and Continued Fraction Reduced Model');

% Frequency Response 
figure;
bode(sys,sys_con);
title('Frequency Response Original and Continued Fraction Reduced Model');
