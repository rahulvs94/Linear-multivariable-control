%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Multivariable Control
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Find the Matrix Sign Function reduced order model, including D*, where gamma is the mean
% of the eigenvalues.
% (b) What are sigma(L1)?
% (c) What are sigma(L2)?
% (d) Plot the output error between the full-order system and the reduced order model given
% step inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
% Given system
A = [0.5 -3  10.5 -21.5  39;
     0.5 -2 -3.5   6.5  -5;
    -1    2 -9     17   -30;
    -0.5  1 -2.5   4.5  -25;
     0    0  0     0    -8];
B = [-2  1  0;
      1 -2 -3;
      1  1  3;
      2  1  4;
      1  0  1];
C = [10 -15 41 -65 131;
     5  -2  10 -2  6];

% Eigenvalues of matrix A
eig(A)

% State space representation
sys = ss(A,B,C,0)

% Shifting the imaginary axis
gamma = mean(eig(A))
A_hat = A-(gamma*eye(5))
eig(A_hat)

% Performing the iterations
sign_lambda = 0.5*(A_hat+inv(A_hat));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda));
sign_lambda = 0.5*(sign_lambda+inv(sign_lambda))


% Getting sign_plus and sign_minus matrices
sign_plus = 0.5*(sign_lambda+eye(5))
sign_minus = 0.5*(eye(5)-sign_lambda)

fprintf('Number of unstable eigenvalues: %d\n',trace(sign_plus));
fprintf('Number of stable eigenvalues: %d\n',trace(sign_minus));

% Getting the matrix of independent column vectors from sign_plus & sign_minus
m1 = sign_plus(:,[1 3 4])
m2 = sign_minus(:,[3 5])

M = [m1 m2]
inv(M)

% Getting matrix sign function
Amsf = inv(M)*A*M
Bmsf =inv(M)*B
Cmsf = C*M

% Lambda 1
lambda_1 = Amsf(1:3,1:3)
eig_lambda1 = eig(lambda_1)
% Lambda 2
lambda_2 = Amsf(4:5,4:5)
eig_lambda2 = eig(lambda_2)

% Reduced order model
Ar = Amsf(1:3,1:3)
Br = Bmsf(1:3,1:3)
Cr = Cmsf(1:2,1:3)
D_star = -C*inv(A)*B + Cmsf*inv(Amsf)*Bmsf

% State space representation of reduced order model
sys_r = ss(Ar,Br,Cr,D_star)

% Output response of original system & reduced order model
figure;
step(sys,'r');
figure;
step(sys_r,'b');
% Error plot
figure;
step(sys-sys_r,'y');
