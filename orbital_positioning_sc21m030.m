clear
clc
%% === Orbital Positioning: Function of Time and Universal Variables === %%
%
% The code presents a method to calculate the orbital position as function
% of time utilizing iterative method complemented by universal variable
% conditionality.
%
% R0        = Initial Position Vector (km)
% V0        = Initial Velocity Vector (km/s)
% R         = Final Position Vector (km)
% V         = Final Velocity Vector (km/s)
% R0m       = Magnitude of Initial Position Vector (km)
% V0m       = Magnitude of Initial Velocity Vector (km/s)
% Rm        = Magnitude of Final Position Vector (km)
% Vm        = Magnitude of Final Velocity Vector (km/s)
% vr0       = Radial Velocity (km/s)
% dt        = Orbital Time Step (sec)
% mu        = Gravitational Parameter (km^3/sec^2)
% aplha     = Reciprocal of semi-major axis (km^-1)
% vr0       = Radial Velocity (km/s)
% xi        = Initial Estimate of Universal Anamoly (km^0.5)
% xf        = Final Estimate of Universal Anamoly (km^0.5)
% c         = C(z) class of Stumpff Functions
% s         = S(z) class of Stumpff Functions
% f         = Lagrange Coefficient f
% g         = Lagrange Coefficient g (sec^-1)
% t         = Trajectory Estimate
% f_xf      = Iterative Equation for universal anamoly estimation
% df_xf     = Iterative Equation change with respect to universal anamoly
% fdot      = Lagrange Coefficient f change with respect to time (sec^-1)
% gfot      = Lagrange Coefficient g change with respect to time
%
% Note:
% 1. The vectors are plotted and considered over i-j-k plane.
% 2. Universal Values are considered in reference to Earth.
% 3. Error Rate for iterative scheme is taken to be 1e-09.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%% ~~~~~~~~~~~~~~~~~~~ Input and Initial Variables ~~~~~~~~~~~~~~~~~~~~~ %%

% Inputing required data for the analysis:

% Input: Position Vector:
R0_i = input('Enter position i-component (|value(i)| > 6378 km): ');
R0_j = input('Enter position j-component (|value(j)| > 6378 km): ');
R0_k = input('Enter position k-component (|value(k)| > 6378 km): ');

% Input: Velocity Vector:
V0_i = input('Enter velocity i-component (km/s): ');
V0_j = input('Enter velocity j-component (km/s): ');
V0_k = input('Enter velocity k-component (km/s): ');

%Input: Time Step:
dt = input('Enter the time step (sec): ');

%Declaring Universal variables and constants:
mu = 398600;

% Forming and calculating position and velocity vector and magnitude
% respectively:
R0 = [R0_i R0_j R0_k];
V0 = [V0_i V0_j V0_k];
R0m = sqrt((R0_i)^2 + (R0_j)^2 + (R0_k)^2);
V0m = sqrt((V0_i)^2 + (V0_j)^2 + (V0_k)^2);

% Calculating for reciprocal of semi-major axis and radial velocity:
alpha = (2/R0m) - ((V0m^2)/mu);
vr0 = (dot(R0, V0))/R0m;

% Calculating for Universal Anamoly:
xi = sqrt(mu)*abs(alpha)*dt;
xf = xi;

% Defining arbitrary constant:
zi = alpha*(xi)^2;

%% ~~~~~~~~~~~~~ Itervative Estimation: Universal Variables ~~~~~~~~~~~~ %%

% Defining required error and ratio for the iterative scheme:
error = 1.e-8;
ratio = 1;

% Defining C(z), S(z) stumpff functions and iterative scheme:
while abs(ratio) > 1.e-9
    zf = alpha*(xf)^2;
    
    % Defininig S(z) class of stumpff functions:
    if zf > 0
        s = (sqrt(zf) - sin(sqrt(zf)))/(sqrt(zf))^3;
    elseif zf < 0
        s = (sinh(sqrt(-zf)) - sqrt(-zf))/(sqrt(-zf))^3;
    else
        s = 1/6;
    end
    
    %Defining C(z) class of stumpff functions:
    if zf > 0
        c = (1 - cos(sqrt(zf)))/zf;
    elseif zf < 0
        c = (cosh(sqrt(-zf)) - 1)/(-zf);
    else
        c = 1/2;
    end
    
    f_xf = (R0m*vr0/sqrt(mu))*(xf^2)*c + (1 - alpha*R0m)*(xf^(3))*s + R0m*xf - sqrt(mu)*dt;
    df_xf = (R0m*vr0/sqrt(mu))*xf*(1 - alpha*(xf^2)*s)+ (1 - alpha*R0m)*(xf^2)*c + R0m;
    
    ratio = f_xf/df_xf;
    xf = xf - ratio;
end

%% ~~~~~~~~~~~~ Output: Final Position and Velocity Vectors ~~~~~~~~~~~~ %%

% Solving for Lagrange Coefficients:
f = 1- ((xf^2)/R0m)*c;
g = dt- (1/sqrt(mu))*(xf^3)*s;

% Calculating for final position vector:
R = f*R0 + g*V0;
Rm = norm(R);

% Solving for Lagrange Coefficients with change in universal anamoly:
fdot = (sqrt(mu)/(Rm*R0m))*(zf*s - 1)*xf;
gdot = 1 - (xf^2/Rm)*c;

% Calculating for final velocity vector:
V = fdot*R0 + gdot*V0;

%Defining trajectory shape:
    if alpha > 0
        t = 'Trajectory is Elliptic';
    elseif alpha < 0
        t = 'Trajectory is Hyperbolic';
    else
        t = 'Trajectory is Parabolic';
    end

% Representing the initial and final outputs of the position and velocity
% vectors:

fprintf('---------------------------------------------------')
fprintf('\n Initial position vector (km):')
fprintf('\n R0 = (%g, %g, %g)\n', R0(1), R0(2), R0(3))
fprintf('\n Initial velocity vector (km/s):')
fprintf('\n V0 = (%g, %g, %g)\n', V0(1), V0(2), V0(3))
fprintf('\n %s', t)
fprintf('\n\n Elapsed time = %g s\n',dt)
fprintf('\n Final position vector (km):')
fprintf('\n R = (%g, %g, %g)\n', R(1), R(2), R(3))
fprintf('\n Final velocity vector (km/s):')
fprintf('\n V = (%g, %g, %g)', V(1), V(2), V(3))
fprintf('\n-----------------------------------------------\n')