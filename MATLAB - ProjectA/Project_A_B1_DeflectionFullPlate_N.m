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
%%%%%%%     B1 - Deflection and maximum stresses, distributed       %%%%%%%
%%%%%%%     static pressure, orthotropic plate, Navier Method       %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

% Steel parameters
E = 2.1e5;           % [MPa]    Young's modulus         
nu = 0.3;            % [-]      Poisson
sigma_y = 250;       % [MPa]    Yielding stress      % NOTE: different from project
%ro = 7850;          % [kg/m^3] Density


%     Table of content of the file:
%
%     Input data
%     Equations definition
%         Navier method - along x2 Find max deflection (-> max M) 
%         Navier method - along x1 Find max deflection (-> max M) 
%         Compute tension and pressure for yelding, for both stiffeners


%% Input data (plate)

a = 55.8*1000;              % [mm] Length
b = 9.95*1000;              % [mm] Width        
t = 6.5;                    % [mm] Thickness
n_1 = floor(b/670);         % Number of stiffeners in x_1 axis
n_2 = floor(a/2800);        % Number of stiffeners in x_2 axis

%% Equations and variable definition

% Geometrical parameters for the stiffeners in the x_1 and x_2 axis
% note: Stiffener 1 is the one along the x_1 direction
f_1  = 12;  f_2 = 100;  % Width of the stiffener [mm]
hf_1 = 10; hf_2 = 10;   % thickness of the stiffener in the horizontal part [mm]
hw_1 = 12; hw_2 = 10;   % same but for vertical part [mm]
h_1  = 155; h_2 = 695;  % high of the stiffener [mm] (called d in the image)
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
d_1 = b/(n_1+1); d_2 = a/(n_2+1); % distance between stiffeners 
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



%% Navier method - Find max deflection (-> max M) - Along x2

syms x_1 x_2 m n
index = 0;

% define the number of iterations
interval = b/50; % Define "every how many points" you compute the momentum
fraction = 0.1;  % define if you want to copmute it for all plate (1) of f.e. half of it (0.5)
limit = 10;
p = 1;            % For computation speed we declare p as 1 (avoid syms)

% Prepare for data collection
data_w = []; sum_w = 0;
data_M_11 = []; sum_M_11 = 0;

% Prepare equations
P_mn = 4/(a*b) * int(int(p * sin((m * pi * x_1)/a) * sin((n * pi * x_2)/b), x_2, 0, b), x_1, 0, a);
w = (P_mn / (pi^4 * ( D_1 * (m / a)^4 + 2*D_12 * (m * n / (a * b))^2 + D_2 * (n / b)^4))) * sin(m * pi * x_1 / a) * sin(n * pi * x_2 / b);
% w_max= vpa(subs(w, [x_1 x_2], [x1_val x2_val]));
M_11 = D_1 * (diff(w,x_1,2) + nu*diff(w,x_2,2));
% M_11_max = vpa(subs(M_11, [x_1 x_2], [x1_val x2_val]));

% Define points for deflection 
x1_val = a/2;
x2_val = 0:interval:b*fraction; % Define vector of points where to calculate deflection 
                  
% Compute deflection in each point
for k=1:length(x2_val)
    %w_max = vpa(subs(w, [x_1 x_2], [x1_val x2_val(k)]));
    M_11_max = vpa(subs(M_11, [x_1 x_2], [x1_val x2_val(k)]));
        for i = 1:limit
        for j = 1:limit
            %sum_w = sum_w + vpa(subs(w_max, [m n], [i j]));
            sum_M_11 = sum_M_11 + vpa(subs(M_11_max, [m n], [i j]));
            %deflection = deflection + vpa(subs(w_max, [m n], [i j]));
        end  
        end   
    index = index+1;
    %data_w(index) = deflection;
    data_M_11(index) = -sum_M_11; % Collect each iteration in an array (- is to fix the sign)
    sum_M_11 = 0;                 % Reset the summation
end

%figure
%plot(1:length(data_w), data_w)
figure
plot(1:length(data_M_11), data_M_11)

% Find peak (maximum deflection) and its location
[M_max,M_max_locs] = findpeaks(data_M_11);
syms p
M_max = M_max*p;                % Make the momentum function of p
M_max_locs = M_max_locs * interval; % Trace back the position in the plate
findpeaks(data_M_11)



%% Navier method - Find max deflection (-> max M) - Along x1

syms x_1 x_2 m n
index = 0;

% define the number of iterations
interval = a/300; % Define "every how many points" you compute the momentum
fraction = 1/10;  % define if you want to copmute it for all plate (1) of f.e. half of it (0.5)
limit = 25;
p = 1;            % For computation speed we declare p as 1 (avoid syms)

% Prepare for data collection
data_w = []; sum_w = 0;
data_M_11 = []; sum_M_11 = 0;

% Prepare equations
P_mn = 4/(a*b) * int(int(p * sin((m * pi * x_1)/a) * sin((n * pi * x_2)/b), x_2, 0, b), x_1, 0, a);
w = (P_mn / (pi^4 * ( D_1 * (m / a)^4 + 2*D_12 * (m * n / (a * b))^2 + D_2 * (n / b)^4))) * sin(m * pi * x_1 / a) * sin(n * pi * x_2 / b);
% w_max= vpa(subs(w, [x_1 x_2], [x1_val x2_val]));
M_11 = D_1 * (diff(w,x_1,2) + nu*diff(w,x_2,2));
% M_11_max = vpa(subs(M_11, [x_1 x_2], [x1_val x2_val]));

% Define points for deflection 
x1_val = 0:interval:a*fraction; % Define vector of points where to calculate deflection 
x2_val = b/2;
                  
% Compute deflection in each point
for k=1:length(x1_val)
    %w_max = vpa(subs(w, [x_1 x_2], [x1_val x2_val(k)]));
    M_11_max = vpa(subs(M_11, [x_1 x_2], [x1_val(k) x2_val]));
        for i = 1:2:limit
        for j = 1:2:limit
            %sum_w = sum_w + vpa(subs(w_max, [m n], [i j]));
            sum_M_11 = sum_M_11 + vpa(subs(M_11_max, [m n], [i j]));
            %deflection = deflection + vpa(subs(w_max, [m n], [i j]));
        end  
        end   
    index = index+1;
    %data_w(index) = deflection;
    data_M_11(index) = -sum_M_11; % Collect each iteration in an array (- is to fix the sign)
    sum_M_11 = 0;                 % Reset the summation
end

%figure
%plot(1:length(data_w), data_w)
figure
plot(1:length(data_M_11), data_M_11)

% Find peak (maximum deflection) and its location
[peaks,locs] = findpeaks(data_M_11);
%[M_max,M_max_locs] = max(peaks);
M_max = max(peaks);
M_max_locs = locs(peaks==M_max);
syms p
M_max = M_max*p;                    % Make the momentum function of p
M_max_locs = M_max_locs * interval; % Trace back the position in the plate
findpeaks(data_M_11)


%% Compute tension and pressure for yelding, for both stiffeners

% Compute the maximum p for yelding - x1 (logitudinals)
b_e_x1 = d_1*min(1.1/(1+2*(d_1/d_2)^2),1);     % effective breadth NOTE: CL from the formula always needs extra care
A_x1 = Astiff_1 + b_e_x1*t;                    % Area of considered stiffener
e_e_x1 = ((-b_e_x1*t*(t/2))+(Astiff_1*e_1))/A_x1;       

I_barred_x1 = (f_1*hf_1^3)/12 + hf_1*f_1*(h_1-e_e_x1)^2 + (hw_1*hb_1^3)/12 + ...
    hw_1*hb_1*((hb_1/2)-e_e_x1)^2 + (b_e_x1*t^3)/12 + b_e_x1*t*((t/2)+e_e_x1)^2;

W_e_x1 = I_barred_x1/((h_1+(hf_1/2)-e_e_x1)*b_e_x1);
sigma_max_x1 = M_max / W_e_x1;
p_x1 = double(solve(sigma_max_x1==sigma_y,p))       % Result in MPa

% Compute the maximum p for yelding - x2 (transversal)
b_e_x2 = d_2*min(1.1/(1+2*(d_2/b)^2),1);            % effective breadth NOTE: CL from the formula always needs extra care
A_x2 = Astiff_2 + b_e_x2*t;                         % Area of considered stiffener
e_e_x2 = ((-b_e_x2*t*(t/2))+(Astiff_2*e_2))/A_x2;   % Same story 

I_barred_x2 = (f_2*hf_2^3)/12 + hf_2*f_2*(h_2-e_e_x2)^2 + (hw_2*hb_2^3)/12 + ...
    hw_2*hb_2*((hb_2/2)-e_e_x2)^2 + (b_e_x2*t^3)/12 + b_e_x2*t*((t/2)+e_e_x2)^2;

W_e_x2 = I_barred_x2/((h_2+(hf_2/2)-e_e_x2)*b_e_x2);
sigma_max_x2 = M_max / W_e_x2;
p_x2 = double(solve(sigma_max_x2==sigma_y,p))       % Result in MPa


