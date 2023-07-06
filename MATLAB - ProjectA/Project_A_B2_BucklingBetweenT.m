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
%%%%%%%   B2 - Buckling of plate in between transveral stiffners    %%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

% Steel parameters
E = 2.1e5;           % [MPa]    Young's modulus         
nu = 0.3;            % [-]      Poisson
sigma_y = 250;       % [MPa]    Yielding stress      
%ro = 7850;          % [kg/m^3] Density


%     Table of content of the file:
%
%     Input data
%     Equations definition
%         Buckling - energy method
%         Find m and sigma_c
%         Johnson - Ostenfeld parabola - Find the corrected sigma_c



%% Input data (plate)

a = 2790;              % [mm] Length
b = 9950;              % [mm] Width        
t = 6.5;                    % [mm] Thickness
n_1 = floor(b/670);         % Number of stiffeners in x_1 axis
n_2 = 0;        % Number of stiffeners in x_2 axis

%% Equations and variable definition

% Geometrical parameters for the stiffeners in the x_1 and x_2 axis
% note: Stiffener 1 is the one along the x_1 direction
f_1  = 12;  f_2 = 0;  % Width of the stiffener [mm]
hf_1 = 10; hf_2 = 0;   % thickness of the stiffener in the horizontal part [mm]
hw_1 = 12; hw_2 = 0;   % same but for vertical part [mm]
h_1  = 155; h_2 = 0;  % high of the stiffener [mm] (called d in the image)
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
K_1 = (1/3) * (hw_1^3*h_1 + hf_1^3*f_1); % In our case, open profile, K is close to zero
K_2 = (1/3) * (hw_2^3*h_2 + hf_2^3*f_2);
G = E/(2*(1+nu));

a_2=0;
% Rigidities for Orthotropic case
D_1 = D * (1 + b_1 + (12 * (1-nu^2) * a_1 * (1+a_2) * c_1^2) / ...
    ((1+a_1) * (1+a_2) - nu^2*a_1*a_2));
D_2 = D * (1 + b_2 + (12 * (1-nu^2) * a_2 * (1+a_1) * c_2^2) / ...
    ((1+a_1) * (1+a_2) - nu^2*a_1*a_2)); 
D_12 = D * (1 + (12*(1-nu^2)*nu*a_1*a_2*c_1*c_2) / ...
    ((1+a_1) * (1+a_2) - nu^2*a_1*a_2)) + G/2*(K_1/d_1 + K_2/d_2);

%D_12 = D;      % Need to impose them, as I am not considering anymore the stiffeners!
%D_2 = D;


%% Buckling - energy method
% longitudinals are smeared out, by copmuting the equivalent thickness h_
% the transversal are instead accounted directly as U_stiff

syms sigma_c x_1 x_2

% Choose how many m and n and built the Wmn consequentially
n = 1;                                % Choose how many n´s!
syms m
% m_val = round(a/b);                   % automatic how many m´s
% if m_val==0  m_val=1;  end            % if closest integer is 0, fix
% m = m_val;

% Compute I_barred_x2 and releated things m
% b_e_x2 = d_2*min(1.1/(1+2*(d_2/b)^2),1);     % effective breadth NOTE: CL from the formula always needs extra care
% A_x2 = Astiff_2 + b_e_x2*t;                  % Area of considered stiffener
% e_e_x2 = ((-b_e_x2*t*(t/2))+(Astiff_2*e_2))/A_x2;        
% I_barred_x2 = (f_2*hf_2^3)/12 + hf_2*f_2*(h_2-e_e_x2)^2 + (hw_2*hb_2^3)/12 + ...
%     hw_2*hb_2*((hb_2/2)-e_e_x2)^2 + (b_e_x2*t^3)/12 + b_e_x2*t*((t/2)+e_e_x2)^2;
% I_f2 = (1/12) * hf_2*f_2^3 + (1/12) * h_2*hw_2^3;

% along x1
% t_ = t + Astiff_1/d_1;   % Compute equivalent thickness

% Load functions
f_11 = t;              % If it is othotropic you need t_ not t !
f_22 = 0;
f_12 = 0;

X_m1 = sin(m*pi*x_1/a); % Simply supported case
X_n2 = sin(n*pi*x_2/b);

w = X_m1*X_n2;

% Compute energies
V = - sigma_c/2 * int(int(f_11*(diff(w,x_1))^2+2*f_12*diff(w,x_1)*...
       diff(w,x_2)+f_22*(diff(w,x_2))^2,x_1,0,a),x_2,0,b);

U =  (1/2) * int(int(D_1*(diff(w,x_1,2))^2 + D_2*(diff(w,x_2,2))^2 + ...
       2*D_12*(diff(w,x_1,2)*diff(w,x_2,2)),x_1,0,a),x_2,0,b);
   
% U_stiff_equ = (1/2) * int(E*I_barred_x2*(diff(w,x_2,2))^2 + G*K_2*(diff(w,x_2)*diff(w,x_1))^2 + E*I_f2*(diff(hb_2*diff(w,x_1),x_2,2))^2,x_2,0,b);
% U_stiff = 0;
% for k=1:n_2
%     U_stiff = U_stiff + subs(U_stiff_equ,x_1,d_2*k);
% end  

% Equilibrium
PI = V + U;



%% Find m and sigma_c

% Find m for minimization
sigma_c = solve(PI==0,sigma_c);        % Create equation for sigma_c from equilibrium
sigma_c = matlabFunction(sigma_c);     % Convert it to function
estimate = a/b;                        % Estimate of m, to find the right area
m_min = fminsearch(sigma_c,estimate);  % Compute local minimum of the function
m = round(m_min)                      % Find the closest integer
if m==0 m = 1; end                        % Condtion if m is 0

sigma_c = double(subs(sigma_c,m))      % Substitute m -> find minimum sigma_c

%% Johnson - Ostenfeld parabola - Find the corrected sigma_c

if sigma_c > 0.5*sigma_y
sigma_c = sigma_y*(1 - sigma_y/(4*sigma_c))
end


