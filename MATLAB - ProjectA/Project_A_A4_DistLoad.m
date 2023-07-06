%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%  41517 Stiffened Plates and Sandwich Constructions  %%%%%%%%%%%
%%%%%%%%%%%%%%%   DTU - TECHNICAL UNIVERSITY OF DENMARK    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Project A: Analysis of a stiffened plate panel  %%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%   Copenhagen, Spring semester 2023   %%%%%%%%%%%
%                                                                         %
%                               Christian Casarotto - s223302             %
%                                    Irene Berganzo - s223230             %
%                                                                         %
%%%%%%%     A4 - deflection and maximum stresses, uniform static    %%%%%%%
%%%%%%%     pressure, isotropic plate                               %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

% Steel parameters
E = 2.1e5;           % [MPa]    Young's modulus         
nu = 0.3;            % [-]      Poisson
sigma_y = 250;       % [MPa]    Yielding stress
ro = 7850;           % [kg/m^3] Density


%     Table of content of the file:
%
%     Input data
%     Equations definition
%     Navier method
%       3D graph (with Navier)
%     Rayleigh-Ritz method
%       Rayleigh-Ritz method for Simply Supported
%       Rayleigh-Ritz method for Clamped B.C.
%       Rayleigh-Ritz method for ...


%% Input data (plate)

a = 2800;        % [mm] Length
b = 670;         % [mm] Width        
t = 6.5;         % [mm] Thickness
n_1 = 0;         % Number of stiffeners in x_1 axis
n_2 = 0;         % Number of stiffeners in x_2 axis

% Load (distributed)
p = -0.01; % MPa 

%% Equations and variable definition

% Geometrical parameters for the stiffeners in the x_1 and x_2 axis
% note: Stiffener 1 is the one along the x_1 direction
f_1 = 0;  f_2 = 0;  % Width of the stiffener [mm]
hf_1 = 0; hf_2 = 0; % thickness of the stiffener in the horizontal part [mm]
hw_1 = 0; hw_2 = 0; % same but for vertical part [mm]
h_1 = 0;  h_2 = 0;  % high of the stiffener [mm] (called d in the image)
hb_1 = h_1-(hf_1/2); hb_2 = h_2-(hf_2/2);

Astiff_1 = f_1*hf_1 + hb_1*hw_1; Astiff_2 = f_2*hf_2 + hb_2*hw_2; % Area of the stiffener [mm^2]
e_1 = (f_1*hf_1*h_1+hb_1^2*hw_1/2) / Astiff_1; e_2 = (f_2*hf_2*h_2+hb_2^2*hw_2/2) / Astiff_2; 
% distance from plate surface to stiffener neutral axis

% Check for NaN
if isnan(e_1) e_1 = 0; end
if isnan(e_2) e_2 = 0; end

I_1 = hw_1*hb_1*(e_1-(hb_1/2))^2 + (hw_1*hb_1^3)/12 + ... 
    hf_1*f_1*(hb_1+(hf_1/2)-e_1)^2 + ((hf_1^3)*f_1)/12;
I_2 = hw_2*hb_2*(e_2-(hb_2/2))^2 + (hw_2*hb_2^3)/12 + ...
    hf_2*f_2*(hb_2+(hf_2/2)-e_2)^2 + ((hf_2^3)*f_2)/12;

% Geometrical parameters of the stiffened plate
D = (E*t^3) / (12*(1-nu^2)); % Isotropic case

% In this case we only have longitudinal stiffeners (along x_1, direction 
% of the ship) therefore all the ones referred to x_2 will be zero (a_2, 
% b_2 ecc) bringing the D2 formulas to just D (automatically adjust).
d_1 = b/(n_1+1); d_2 = a/(n_2+1); % distance between stiffeners along x_1
c_1 = (e_1/t)+0.5; c_2 = (e_2/t)+0.5; % Eccentricity ratio
b_1 = (E*I_1)/(D*d_1); b_2 = (E*I_2)/(D*d_2); % rigidity ratio
a_1 = Astiff_1/(d_1*t); a_2 = Astiff_2/(d_2*t); % are ratio (stiffener/plate)
K_1 = 0; K_2 = 0; % In our case, open profile, K is close to zero
G = E/(2*(1+nu));

% Rigidities for Orthotropic case
D_1 = D * (1 + b_1 + (12 * (1-nu^2) * a_1 * (1+a_2) * c_1^2) / ...
    ((1+a_1) * (1+a_2) - nu^2*a_1*a_2));
D_2 = D * (1 + b_2 + (12 * (1-nu^2) * a_2 * (1+a_1) * c_2^2) / ...
    ((1+a_1) * (1+a_2) - nu^2*a_1*a_2)); 
D_12 = D * (1 + (12*(1-nu^2)*nu*a_1*a_2*c_1*c_2) / ...
    ((1+a_1) * (1+a_2) - nu^2*a_1*a_2)) + G/2*(K_1/d_1 + K_2/d_2);

%% Navier method

syms x_1 x_2 
w_tot = 0; data_w = []; index = 0;

% define the number of iterations
limit = 20;

syms m n

sum_w = 0;

% x1_val = a/2;
% x2_val = b/2;
P_mn = 4/(a*b) * int(int(p * sin((m * pi * x_1)/a) * sin((n * pi * x_2)/b), x_2, 0, b), x_1, 0, a);
w = (P_mn ./ (pi.^4 .* ( D_1 .* (m ./ a).^4 + 2.*D_12 .* (m .* n ./ (a .* b)).^2 + D_2 .* (n ./ b).^4)))  .* sin(m .* pi .* x_1 ./ a) .* sin(n .* pi .* x_2 ./ b);
w_max= vpa(subs(w, [x_1 x_2], [a/2 b/2]));

% Compute stresses
x_3 = t/2;
sigma_1 = -(E*x_3/(1-nu^2))*(diff(w,x_1,2)+nu*diff(w,x_2,2));
sigma_1_point = vpa(subs(sigma_1, [x_1 x_2], [a/2 b/2]));
sigma_2 = -(E*x_3/(1-nu^2))*(diff(w,x_2,2)+nu*diff(w,x_1,2));
sigma_2_point = vpa(subs(sigma_2, [x_1 x_2], [a/2 b/2]));
tau_12 = -2*G*x_3*diff(w,x_1)*diff(w,x_2);
tau_12_point = vpa(subs(tau_12, [x_1 x_2], [a/2 b/2]));

data_sigma_1 = [];
data_sigma_2 = [];
data_tau_12 = [];
data_sigma_VM = [];
sum_sigma_1 = 0;
sum_sigma_2 = 0;
sum_tau_12 = 0;

% compute double summation
for i = 1:2:limit
    for j = 1:2:limit

      sum_w = sum_w + vpa(subs(w_max, [m n], [i j]));

      sigma_1_values = vpa(subs(sigma_1_point, [m n], [i j]));
      sum_sigma_1 = sum_sigma_1 + sigma_1_values;

      sigma_2_values = vpa(subs(sigma_2_point, [m n], [i j]));
      sum_sigma_2 = sum_sigma_2 + sigma_2_values;

      tau_12_values = vpa(subs(tau_12_point, [m n], [i j]));
      sum_tau_12 = sum_tau_12 + tau_12_values;

      sigma_VM = sqrt(sum_sigma_1^2 + sum_sigma_2^2 -sum_sigma_1*sum_sigma_2 + 3*sum_tau_12);

      % Collect each iteration in an array
      index = index+1;
      data_w(index) = sum_w;
      data_sigma_1(index) = sum_sigma_1;
      data_sigma_2(index) = sum_sigma_2;
      data_tau_12(index) = sum_tau_12;
      data_sigma_VM(index) = sigma_VM;

    end  
end

figure
plot(1:length(data_w), data_w)
figure
plot(1:length(data_sigma_1), data_sigma_1)
figure
plot(1:length(data_sigma_2), data_sigma_2)
figure
plot(1:length(data_tau_12), data_tau_12)
figure
plot(1:length(data_sigma_VM), data_sigma_VM)



%% 3D graph with Navier method

syms x_1 x_2 
syms m n
x_3 = t/2;

% define the number of iterations and number of points per edge
limit = 13;  % the functions converge around 13 iterations
np = 10;     % for the pictures in the report we used 80 points

% Definition of Navier functions
P_mn = 4/(a*b) * int(int(p * sin((m * pi * x_1)/a) * sin((n * pi * x_2)/b), x_2, 0, b), x_1, 0, a);
w = (P_mn / (pi^4 * ( D_1 * (m / a)^4 + 2*D_12 * (m * n / (a * b))^2 + D_2 * (n / b)^4)))  * sin(m * pi * x_1 / a) * sin(n * pi * x_2 / b);

% Compute stresses 
sigma_1 = -(E*x_3/(1-nu^2))*(diff(w,x_1,2)+nu*diff(w,x_2,2));
sigma_2 = -(E*x_3/(1-nu^2))*(diff(w,x_2,2)+nu*diff(w,x_1,2));
tau_12 = -2*G*x_3*diff(w,x_1)*diff(w,x_2);

Cx1 = (linspace(0,a,np));
Cx2 = (linspace(0,b,np))';

% Deflection
Deflection=0;
Deflection1D = subs(w,[x_1],[Cx1]);
Deflection2D = subs(Deflection1D,[x_2],[Cx2]);
for i=1:2:limit
    for j=1:2:limit
        DeflectionC = subs(Deflection2D,[m n],[i j]);
        Deflection = Deflection + DeflectionC;
    end
end
Deflection=double(Deflection);

% Sigma_1
sigma_1_R=0;
sigma_11D = subs(sigma_1,[x_1],[Cx1]);
sigma_12D = subs(sigma_11D,[x_2],[Cx2]);
for i=1:2:limit
    for j=1:2:limit
        sigma_1C = subs(sigma_12D,[m n],[i j]);
        sigma_1_R = sigma_1_R + sigma_1C;
    end
end
sigma_1_R=double(sigma_1_R);

% Sigma_2
sigma_2_R=0;
sigma_21D = subs(sigma_2,[x_1],[Cx1]);
sigma_22D = subs(sigma_21D,[x_2],[Cx2]);
for i=1:2:limit
    for j=1:2:limit
        sigma_2C = subs(sigma_22D,[m n],[i j]);
        sigma_2_R = sigma_2_R + sigma_2C;
    end
end
sigma_2_R=double(sigma_2_R);

% tau_12
tau_12_R=0;
tau_121D = subs(tau_12,[x_1],[Cx1]);
tau_122D = subs(tau_121D,[x_2],[Cx2]);
for i=1:2:limit
    for j=1:2:limit
       tau_12C = subs(tau_122D,[m n],[i j]);
       tau_12_R = tau_12_R + tau_12C;
    end
end
tau_12_R=double(tau_12_R);

% sigma_VM
sigma_VM_R = sqrt(abs(sigma_1_R.^2 + sigma_2_R.^2 -sigma_1_R.*sigma_2_R + 3.*tau_12_R));
sigma_VM_R = double(sigma_VM_R);

% Plotting
colormap(figure,winter);
surf(Cx1,Cx2,Deflection)
title('Deflection')
xlabel('x1 [mm]'), ylabel('x2 [mm]'), zlabel('w(x1,x2) [mm]')
ax = gca; 
ax.FontSize = 14; 
colorbar
colormap(figure,winter);
surf(Cx1,Cx2,sigma_1_R)
title('Stress in x_1 direction')
xlabel('x1 [mm]'), ylabel('x2 [mm]'), zlabel('\sigma_1 [MPa]')
ax = gca; 
ax.FontSize = 14;
colorbar
colormap(figure,winter);
surf(Cx1,Cx2,sigma_2_R)
title('Stress in x_2 direction')
xlabel('x1 [mm]'), ylabel('x2 [mm]'), zlabel('\sigma_2 [MPa]')
ax = gca; 
ax.FontSize = 14;
colorbar
colormap(figure,winter);
surf(Cx1,Cx2,tau_12_R)
title('Shear stress')
xlabel('x1 [mm]'), ylabel('x2 [mm]'), zlabel('\tau_{12} [MPa]')
ax = gca; 
ax.FontSize = 14;
colorbar
colormap(figure,winter);
surf(Cx1,Cx2,sigma_VM_R)
title('Von Mises stress')
xlabel('x1 [mm]'), ylabel('x2 [mm]'), zlabel('\sigma_{VM} [MPa]')
ax = gca; 
ax.FontSize = 14;
colorbar



%% Rayleigh-Ritz method simply supported

syms x_1 x_2 m n
w_tot = 0;
data = [];
index = 0;

% Define points where to compute deflection
x1_val = a/2; 
x2_val = b/2;
sum_w = 0;

% Define the maximum m and n
limit = 20;

% Define shape functions for simply supported and deflection
syms W
X_m1 = sin(m*pi*x_1/a);
X_n2 = sin(n*pi*x_2/b);    
%X_m1 = 1 - cos(m*2*pi*x_1/a);
%X_n2 = 1 - cos(n*2*pi*x_2/b);
w = W * X_m1 * X_n2;

% Energies
V = - int(int(p*w,x_1,0,a),x_2,0,b);
U = (1/2) * int(int(D_1*(diff(w,x_1,2))^2 + D_2*(diff(w,x_2,2))^2 + ...
      2*D_12*(diff(w,x_1,2)*diff(w,x_2,2)),x_1,0,a),x_2,0,b);

% Equilibrium
PI = V + U;

% Solve for coefficient W
W = vpa(solve(diff(PI,W)==0,W));

% Compute deflection
w = W * X_m1 * X_n2;
w_max= vpa(subs(w, [x_1 x_2], [2800 670]));

% Compute stresses
x_3 = t/2;
sigma_1 = -(E*x_3/(1-nu^2))*(diff(w,x_1,2)+nu*diff(w,x_2,2));
sigma_1_point = vpa(subs(sigma_1, [x_1 x_2], [a/2 b/2]));
sigma_2 = -(E*x_3/(1-nu^2))*(diff(w,x_2,2)+nu*diff(w,x_1,2));
sigma_2_point = vpa(subs(sigma_2, [x_1 x_2], [a/2 b/2]));
tau_12 = -2*G*x_3*diff(w,x_1)*diff(w,x_2);
tau_12_point = vpa(subs(tau_12, [x_1 x_2], [2800 670]));

data_sigma_1 = [];
data_sigma_2 = [];
data_tau_12 = [];
data_sigma_VM = [];
sum_sigma_1 = 0;
sum_sigma_2 = 0;
sum_tau_12 = 0;
sum_sigma_VM = 0;

% compute double summation
for i = 1:limit
  for j = 1:limit

      w_trial = vpa(subs(w_max, [m n], [i j]));
      sum_w = sum_w + w_trial;

      sigma_1_values = vpa(subs(sigma_1_point, [m n], [i j]));
      sum_sigma_1 = sum_sigma_1 + sigma_1_values;

      sigma_2_values = vpa(subs(sigma_2_point, [m n], [i j]));
      sum_sigma_2 = sum_sigma_2 + sigma_2_values;

      tau_12_values = vpa(subs(tau_12_point, [m n], [i j]));
      sum_tau_12 = sum_tau_12 + tau_12_values;

      sigma_VM = sqrt(sum_sigma_1^2 + sum_sigma_2^2 -sum_sigma_1*sum_sigma_2 + 3*sum_tau_12);

      % Collect each iteration in an array
      index = index+1;
      data_w(index) = sum_w;
      data_sigma_1(index) = sum_sigma_1;
      data_sigma_2(index) = sum_sigma_2;
      data_tau_12(index) = sum_tau_12;
      data_sigma_VM(index) = sigma_VM;

  end  
end

figure
plot(1:length(data_w), data_w)
figure
plot(1:length(data_sigma_1), data_sigma_1)
figure
plot(1:length(data_sigma_2), data_sigma_2)
figure
plot(1:length(data_tau_12), data_tau_12)
figure
plot(1:length(data_sigma_VM), data_sigma_VM)



%% Rayleigh-Ritz method for Clamped B.C.

syms x_1 x_2 m n
w_tot = 0;
data = [];
index = 0;

% Define points where to compute deflection
x1_val = a/2; 
x2_val = b/2;
sum = 0;

% Define the maximum m and n
limit = 20;

% Define shape functions for simply supported and deflection
syms W
X_m1 = 1 - cos(m*2*pi*x_1/a);
X_n2 = 1 - cos(n*2*pi*x_2/b);
%X_m1 = sin(m*pi*x_1/a);      % Simply supported case
%X_n2 = sin(n*pi*x_2/b);    
w = W * X_m1 * X_n2;

% Energies
V = - int(int(p*w,x_1,0,a),x_2,0,b);
U = (1/2) * int(int(D_1*(diff(w,x_1,2))^2 + D_2*(diff(w,x_2,2))^2 + ...
      2*D_12*(diff(w,x_1,2)*diff(w,x_2,2)),x_1,0,a),x_2,0,b);

% Equilibrium
PI = V + U;

% Solve for coefficient W
W = vpa(solve(diff(PI,W)==0,W));

% Compute deflection
w = W * X_m1 * X_n2;
w_max= vpa(subs(w, [x_1 x_2], [a/2 b/2]));

% compute double summation
for i = 1:limit
  for j = 1:limit

      w_trial = vpa(subs(w_max, [m n], [i j]));
      sum = sum + w_trial;

      % Collect each iteration in an array
      index = index+1;
      data(index) = sum;

  end  
end

figure
plot(1:length(data), data)



%% Rayleigh-Ritz method for Clamped along b and Simply Supported along a

syms x_1 x_2 m n
w_tot = 0;
data = [];
index = 0;

% Define points where to compute deflection
x1_val = a/2; 
x2_val = b/2;
sum = 0;

% Define the maximum m and n
limit = 20;

% Define shape functions for simply supported and deflection
syms W
%X_m1 = 1 - cos(m*2*pi*x_1/a);   % clamped case
X_n2 = 1 - cos(n*2*pi*x_2/b);
X_m1 = sin(m*pi*x_1/a);          % Simply supported case
%X_n2 = sin(n*pi*x_2/b); 
w = W * X_m1 * X_n2;

% Energies
V = - int(int(p*w,x_1,0,a),x_2,0,b);
U = (1/2) * int(int(D_1*(diff(w,x_1,2))^2 + D_2*(diff(w,x_2,2))^2 + ...
      2*D_12*(diff(w,x_1,2)*diff(w,x_2,2)),x_1,0,a),x_2,0,b);

% Equilibrium
PI = V + U;

% Solve for coefficient W
W = vpa(solve(diff(PI,W)==0,W));

% Compute deflection
w = W * X_m1 * X_n2;
w_max= vpa(subs(w, [x_1 x_2], [a/2 b/2]));

% compute double summation
for i = 1:limit
  for j = 1:limit

      w_trial = vpa(subs(w_max, [m n], [i j]));
      sum = sum + w_trial;

      % Collect each iteration in an array
      index = index+1;
      data(index) = sum;

  end  
end

figure
plot(1:length(data), data)
