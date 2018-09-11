%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Multivariable Control
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Determinent(SI-A)
% (b) Adjoint(SI-A)
% (c) Inverse of (SI-A)
% (d) Y(s) as RMFD
% (e) Eigenvalues
% (f) Modal matrix (M)
% (g) The diagonal linear representation
%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
A = [-0.9 6.75 0.2 0;-1.67 -2.1 0 0.01;0 0 -0.2 -1.2;0.02 0.1 3.33 -0.2]   %% Matrix A
B = [2.5 -0.25;-0.56 1.17;1 0;-1.67 -1.67]                                 %% Matrix B
C = [1 1 0 2;0 1 1 0]                                                      %% Matrix C    
syms s;                                                                    %% Directly assigning symbol
P = s*eye(4)-A;                                                             
determinent = det(P)                                                       %% Determinent of (SI-A)
adjoint = adjoint(P)                                                       %% Adjoint of (SI-A)        
inverse = adjoint/determinent                                              %% Inverse of (SI-A)
system = simplify(C*(inverse)*B)                                           %% 2 input 2 output system transfer function
[ss_n,ss_d] = numden(simplify(system));                                    %% Simplified transfer function
n_a = fliplr(coeffs(ss_n(1,1)));
n_b = fliplr(coeffs(ss_n(1,2)));
n_c = 2*fliplr(coeffs(ss_n(2,1)));
n_d = 2*fliplr(coeffs(ss_n(2,2)));
d = fliplr(coeffs(ss_d(1,1)));                                             %% no need to calculate for all as denominator will be same
% Numerator matrix for RMFD
N1 = [n_a(1) n_b(1);n_c(1) n_d(1)]/d(1)
N2 = [n_a(2) n_b(2);n_c(2) n_d(2)]/d(1)
N3 = [n_a(3) n_b(3);n_c(3) n_d(3)]/d(1)
N4 = [n_a(4) n_b(4);n_c(4) n_d(4)]/d(1)
% Denominator matrix for RMFD
D1 = (d(1)*eye(2))/d(1)
D2 = (d(2)*eye(2))/d(1)
D3 = (d(3)*eye(2))/d(1)
D4 = (d(4)*eye(2))/d(1)
D5 = (d(5)*eye(2))/d(1)
[eigenvectors, eigenvalues] = eig(A)                                       %% Eigenvalues and eigenvectors   

% check below line
eigenvalues_dup = roots(fliplr(coeffs(determinent)));

M = eigenvectors                                                           %% Modal matrix       
D = inv(M)*A*M                                                             %% Diagonal linear representation