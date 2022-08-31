%% Initial weight estimates
% Mission profile
% 0 - Start of take off
% 1 - End of takeoff/start of climb
% 2 - End of climb/start of cruise
% 3 - End of cruise/start of desecent
% 4 - End of descent and have landed

% P0_elec_P0 = 0 % Using fully fuel
P0_elec_P0 = 0.10 % Proportion of power produced by electrical motor
P0_ICE_P0 = 1-P0_elec_P0


% Weight fractions (Raymer, assuming purely fuel burning aircraft)
W1_W0 = 0.97 % Note that if we are using electricity for this section then W1_W0 = 1
W2_W1 = 0.985
W4_W3 = 0.995

% Apply electric proportion
W1_W0 = 1-P0_ICE_P0*(1-W1_W0)
W2_W1 = 1-P0_ICE_P0*(1-W2_W1)
W4_W3 = 1-P0_ICE_P0*(1-W4_W3)


% Gravity
g = 9.8065

% Assume using turboprop
A = 12 % Qunhy's assumption for regional turboprops
% Cd0 = 0.0478; K = 1/33.67 % (Rajan 2006)
% Cl_Cd_max = 1/(4*Cd0*K)^(1/2)
% e = 1/K*1/(pi*A) % Oswald efficiency factor (good to check)

Cd0_to = 0.02, e_to = 0.7; % Qunhy's values
Cl_Cd_max = ((pi*A*e_to)/(4*Cd0_to))^(1/2)

% Require a harmonic range of 700nm and V_cruise = 300KTAS
V_cruise = 300 %KTAS
V_cruise = V_cruise/1.94384 %m/s
R = 700 %nm
R = R*1.852*1000 %m

% Altitude
h = 25000 %ft
h = h*0.3048 %m

% Speed of sound
a_a0 = 0.910112
a0 = 340.3 %m/s
a = a_a0*a0 %m/s

% Mach number
M_cruise = V_cruise/a

% Heating value
H = 42.5 %MJ/kg % Jet A1
H = H*10^6
% Total efficiency
n_total = 0.32 % Set 12 Slide 4
% Propulsive efficiency
n_pr_to = 0.85
% cp = n_pr_to/n_total*1/H
cp = 0.085*10^-6
cp = 0.7*cp

% Calculating W2_W3
W3_W2 = exp(-R*g*cp/(n_pr_to*Cl_Cd_max))
% Calculating weight fraction from 0 to 4
W4_W0 = W1_W0*W2_W1*W3_W2*W4_W3
% Fuel weight fraction
Wf_W0 = (1-W4_W0)
% Empty weight mass fraction
A2 = 0.92 % For turboprop (not aspect ratio) Set 2 slide 4
C = -0.05  % For turboprop Set 2 slide 4

% Initial guess of W0
W0 = 20*10^3 %kg

% Payload weight
Wp = 6300 %kg
% Wp = (82 + (18+27)/2)*54

iterate = 5;
% Iterate

for i = 1:iterate
    % Calculate empty weight mass fraction
    We_W0 = A2*W0^C
    W0 = Wp/(1-Wf_W0-We_W0)
    % Calculate the fuel weight
    Wf = Wf_W0*W0
    % Add the reserve fuel
    W0 = 0.05*Wf + W0
end

%% Mission profile calculation version 1 (here I assumed that the entire takeoff is governed by electric - not feasible)
clear all; close all; clc
% Mission profile
% 0 - Start of take off
% 1 - End of takeoff/start of climb
% 2 - End of climb/start of cruise
% 3 - End of cruise/start of desecent
% 4 - End of descent and have landed

% Takeoff will be powered electrically for the run but when its in the
% air, it will be using the fuel
% All other segments of the flight will use fuel and will recharge the
% batteries

% s_to = s_run + s_air and s_to = 1000m

% Finding s_run (does not use correlations for TOFL)

% m*a = T-D-F_r
% Applying reduction factor to eliminate D and F_r
% m*a_eff = T where a_eff = a/r_T, r_T varies between 0.8 to 0.9 for a jet
% aircraft
% Break up thrust, T, into two components
    % T = (2*P_s^2*rho*A_p)^(1/3) from V = 0 to V_intersect (this is to
    % account for the fact that T is infinite for V = 0
    % T = n_pr*n_e*P_s/V for V between 10 to V2 (for electric)
% a_eff = T/m
% Can obtain an average acceleration through averaging a_eff
    % a_av = int(a_eff, V = 0, V = V2)/(V2-0)
% V2 is given by (2/rho*W/S*1/Cl_2)^(1/2) where Cl_2 = Cl_max/k_to^2 and k_to = 1.2
% Cl_max_to for regional turboprops lie b/w 1.7 to 2.1
% Also P_s = P_s_battery P_s_battery = max_power_output*m_batteries


rho_SL = 1.225 %kg/m^3
W0 = 21000 %kg
W0_S_to = 3900 % Pa (from set 8 slide 2)
Cl_max_to = 1.9 % Averaged
k_to = 1.2
% Calculating V2
Cl_2 = Cl_max_to/k_to^2
V_2 = (2/rho_SL*W0_S_to*1/Cl_2)^(1/2) %m/s
% Calculating propeller diameter
K_p = 0.56 % Assuming two blades (Set 12 Slide 10)

syms m_battery V
assume(m_battery,'positive')
% Electric motor power output
max_power_output = 1.33*10^3 % W/kg
% Calculating first part of thrust
P_s_0_battery = max_power_output*m_battery
D = K_p*P_s_0_battery^(1/4)
A_p = (pi*D^2/4)
T_seg_1 = (2*P_s_0_battery^2*rho_SL*A_p)^(1/3)

% Calculating second part of thrust
n_pr_to = 0.93
n_e = 0.95 % Wk3 Extra Content Set Slide 9
T_seg_2 = n_pr_to*n_e*P_s_0_battery/V

% Finding V_intersect
V_intersect = solve(T_seg_2 == T_seg_1,V)

% Integrate the thrust over V = 0 to V = V2

T_integrated = int(T_seg_1,V,[0,V_intersect]) + int(T_seg_2,V,[V_intersect,V_2])
a_average = T_integrated/W0*1/(V_2)
% Apply reduction factor
r_t = 0.8
a_average = r_t*a_average 

s_run = V_2^2/(2*a_average)


% Attempt 1
% Finding s_air
% To find s_air we need to find gamma_2. To do this we will assume we are
% flying at a maximum rate of climb
% v_max = 1/(L/D)_max*(27/16)^(1/4)*(P_A_av/P_min-1) where v_max =
% (Vv/V_star)_max and V_star = (2/rho*W/S)^(1/2)*(1/Cl_star)^(1/2)
% P_min = (256/27)^(1/4)(2/rho*W/S)^(1/2)*(Cd0K^3)*W
% P_A = n_pr*alpha*P_s_0, P_s_0 = P0, alpha = (sigma-c)/(1-c) and c = 0.132
% sigma = (T/T0)^(-g/(a*R)-1), a = -6.5*10^-3, T0 = 288.16K and T = a*h+T0
% where h is altitude in m
% P_A_av = 1/(h_cruise-0)*int(n_pr*alpha*P_s_0, h = 0, h = h_cruise) and
% P0_W0 = 0.33 % For twin turboprop (Set 12, Slide 2)
% Finally sin(gamma_2) = Vv_max/V where V = V_2

% Attempt 2
% Finding s_air
% s_air is given by h_to/tan(gamma_2)
% sin(gamma_2) = (T/W-Cd/Cl)_2
% We can find T_2 by evaluating T_seg_2 at V2, W_2 = W0, Cl_2 =
% Cl_max_to/kto^2 and C_d = Cd0 + K*Cl_2^2

T_2 = subs(T_seg_2,V,V_2)
% Aerodynamic parameters
Cd0_to = 0.02, e_to = 0.7; % Qunhy's values
A = 12 % Qunhy's assumption for regional turboprops
K_to = 1/(pi*A*e_to)

Cd_2 = Cd0_to + K_to*Cl_2^2

g = 9.8065
sin_gamma_2 = (T_2/(W0*g)-Cd_2/Cl_2)

gamma_2 = asin(sin_gamma_2)

h_to = 35 %ft
h_to = h_to*0.3048 %m

s_air = h_to/tan(gamma_2)

s_to = s_run + s_air
assume(s_run,'positive')
assume(s_air,'positive')


% Plot s_run against s_air
figure
fplot(m_battery,s_run,[0,2000])
hold on
fplot(m_battery,s_air,[0,2000])
hold on
fplot(m_battery,s_to,[0,2000])
hold on
yline(1000)
m_battery_array = linspace(1000,2000,10000);
s_to_array = double(subs(s_to,m_battery,m_battery_array));

% Require battery
dist_array = abs(s_to_array - 1000);
[min_dist,min_index] = min(dist_array,[],'all');
m_battery_estimate = m_battery_array(min_index)


% m_battery_estimate = solve(s_to == 1000,m_battery)
% 
% Checking sensibility of values
gamma_estimate = double(subs(gamma_2,m_battery,m_battery_estimate))*180/pi
V_intersect_estimate = double(subs(V_intersect,m_battery,m_battery_estimate))
figure
fplot(V,subs(T_seg_2,m_battery,m_battery_estimate),[double(V_intersect_estimate),V_2])
hold on
fplot(V,subs(T_seg_1,m_battery,m_battery_estimate),[0,double(V_intersect_estimate)])
a_average_estimate = double(subs(a_average,m_battery,m_battery_estimate))
V_2_estimate = double(subs(V_2,m_battery,m_battery_estimate))
T_2_estimate = double(subs(T_2,m_battery,m_battery_estimate))
T_2_W0 = double(subs(T_2/(W0*g),m_battery,m_battery_estimate))
Cd_Cl_2 = Cd_2/Cl_2
P_s_0_battery_estimate = double(subs(P_s_0_battery,m_battery,m_battery_estimate))
s_run_estimate = double(subs(s_run,m_battery,m_battery_estimate))
s_air_estimate = double(subs(s_air,m_battery,m_battery_estimate))

% Calculate the time taken
time_run = (2*s_run_estimate/a_average_estimate)^(1/2)
time_air = h_to/(sin(gamma_estimate*pi/180)*V_2_estimate) % Assume thrust is relatively constant over air time
time_TO = time_run + time_air

energy_batteries = P_s_0_battery_estimate*m_battery_estimate


%% Climb (maximizing climb rate) (assumed that the climb rate is not the same, not a conventional way of doing things)
syms h
assume(h,'positive')

W1_W0 = 1
W2_W1 = 0.985

T = -6.5*10^-3*h + 288.16 %K
sigma = (T/288.16)^(-g/(-6.5*10^-3*287)-1)
alpha = (sigma-0.132)/(1-0.132)
P0_W0_tp = 0.33 % W/g
P0_tp = P0_W0_tp*W0*10^3 %W
P_s_0_tp_climb = P0_tp

S_climb = W0_S_to^-1*W0*g

% Total efficiency
n_total_climb = 0.32 % Set 12 Slide 4
n_pr_climb = 0.93;

H = 42.5*10^6 %J/g

cp_climb = n_pr_climb/n_total_climb*1/H

Cd0_climb = 0.02; K_climb = K_to;


P_A_climb = n_pr_climb*alpha*P_s_0_tp_climb
% W_climb = W0-cp_climb*alpha*P_s_0_tp_climb*h/Vv_max <- Will need to solve
% an integrodifferetial equation
W_climb = (W0 + (W2_W1)*W1_W0*W0)/2 % Take the average between the start of climb and end of climb as the weight
P_min_climb = (256/27)^(1/4)*(2/(rho_SL*sigma)*W_climb*g/S_climb)^(1/2)*(Cd0_climb*K_climb^3)^(1/4)*W_climb*g

Cl_star = (Cd0_climb/K_climb)^(1/2)

V_star = (2/(rho_SL*sigma)*W_climb*g/S_climb)^(1/2)*(1/Cl_star)^(1/2)

L_D_star_climb = 1/(4*Cd0_climb*K_climb)^(1/2)

Vv_max = 1/L_D_star_climb*(27/16)^(1/4)*(P_A_climb/P_min_climb-1)*V_star

% Results over climb
fuel_burned = W0 - W2_W1*W1_W0*W0
% figure
% fplot(Vv_max,[0 10000])
total_time_taken = vpaintegral(1/Vv_max,h,[0 25000*0.3048])

%% Constraint Analysis
% Mission profile
% 0 - Start of take off
% 1 - End of takeoff/start of climb
% 2 - End of climb/start of cruise
% 3 - End of cruise/start of desecent
% 4 - End of descent and have landed

% Figure formatting
set(0,'DefaultLineLineWidth',1,...+2
    'DefaultLineMarkerSize',10,...
    'DefaultAxesFontSize',16,...
    'DefaultTextInterpreter','latex',...
    'DefaultLegendInterpreter','latex',...
    'DefaultLegendFontSize',16)

lower_W0_S = 1000
upper_W0_S = 8000
resolution = 100

figure_position = [100 100 1000 800]
constraint_plot = figure
constraint_plot.Position = figure_position

rho_SL = 1.225

P0_elec_P0 = 0.1 % Using fully fuel
% P0_elec_P0 = 0.10 % Proportion of power produced by electrical motor
P0_ICE_P0 = 1-P0_elec_P0

% Initial parameters
W1_W0 = 0.97 % Note that if we are using electricity for this section then W1_W0 = 1
W2_W1 = 0.985
W4_W3 = 0.995

% Apply electric proportion
W1_W0 = 1-P0_ICE_P0*(1-W1_W0)
W2_W1 = 1-P0_ICE_P0*(1-W2_W1)
W4_W3 = 1-P0_ICE_P0*(1-W4_W3)

% Note We_W0 = 0.92*M0^-0.05

Cd0 = 0.02 % Torenbeek

Cd0_to = Cd0
Cd0_climb = Cd0
Cd0_cruise = Cd0
Cd0_land = Cd0

% A = 12 % Qunhy's assumption
A = 11.5 % Refined valued
e = 0.7 % Set 6 Slide 2

K = 1/(pi*A*e)

Cl_max_clean = 1.6
Cl_max_to = 2.0
Cl_max_land = 2.3 % Torenbeek ?

n_pr = 0.85
HEP_cp_reduction = 0.3 % If hybrid parallel
% HEP_cp_reduction = 0 % If fuel only
cp_cruise = 0.085*10^-6 % kg/Ws
cp_cruise = cp_cruise*(1-HEP_cp_reduction)
cp_loiter = 0.101*10^-6 %kg/Ws
cp_loiter = cp_loiter*(1-HEP_cp_reduction)

g = 9.8065 %m/s^2

L_D_star = 1/(4*Cd0_cruise*K)^(1/2)

% W3_W2 = exp(-R_h*g*cp_cruise/(n_pr*L_D_star)) % Weight fraction if you fly at V_star

% Analytical curve for sigma
syms h
T = -6.5*10^-3*h + 288.16 %K
sigma = (T/288.16)^(-g/(-6.5*10^-3*287)-1)


% Landing approach speed 
% Assume approaching at Vapp = 1.3*V_stall 
% Also assume landing at SL
% Assume beta_land = 1
k_land = 1.3
% V_app = 55 % Set 6 Slide 1
V_app = 77 %m/s https://skybrary.aero/articles/approach-speed-categorisation
sigma_land = 1
beta_land = 1 % Not sure about this assumption

W0_S_land = (V_app/k_land)^2*(sigma_land*rho_SL/2)*1/beta_land*Cl_max_land % Pa

xline(W0_S_land)
hold on

% Take off field length
% Assume taking off at SL
% TOFL = 1000m (for short take off)
syms W0_S P0_W0
assume(W0_S,'positive')
assume(P0_W0,'positive')

sigma_to = 1

TOP = W0_S*(P0_W0)^-1*1/sigma_to*1/Cl_max_to
TOFL = 11.8*TOP+0.255*TOP^2 % For two engine prop
TOFL_req = 1000

W0_S_to = linspace(lower_W0_S,upper_W0_S,resolution);
P0_W0_to = []
for iW0_S = 1:length(W0_S_to)
    P0_W0_to_temp = solve(TOFL_req == subs(TOFL,W0_S,W0_S_to(iW0_S)),P0_W0);
    P0_W0_to(iW0_S) = double(P0_W0_to_temp);
end

plot(W0_S_to,P0_W0_to)

% Cruise speed

% Alitude density effect
syms h W3_W2
T = -6.5*10^-3*h + 288.16 %K
alpha_cruise = 0.7*sigma^0.8

% Required harmonic range
R_h = 700 %nm
R_h = R_h*1.852*1000 %m

% Cruise speed required
V_cruise_req = 300 %KTAS
V_cruise_req = V_cruise_req/1.94384

% Zero drag
Cd0_cruise = 0.02
% Max lift to drag ratio
L_D_star_cruise = 1/(4*Cd0_cruise*K)^(1/2)

% Cruise weight fraction
beta_cruise = W1_W0*W2_W1*(W3_W2+1)/2

% Lift to drag ratio
Cl_cruise = 2/(sigma*rho_SL*V_cruise_req^2)*beta_cruise*W0_S
Cd_cruise = Cd0_cruise + K*Cl_cruise^2
L_D_cruise = Cl_cruise/Cd_cruise

% Gravity
g = 9.8065

eq0 = W3_W2 == exp(-R_h*g*cp_cruise/(n_pr*L_D_cruise))
% Max lift to drag speed
Cl_star = (Cd0_cruise/K)^(1/2)
V_star = (2/(sigma*rho_SL)*beta_cruise*W0_S)^(1/2)*(1/Cl_star)^(1/2)
u = V_cruise_req/V_star
P0_W0_cruise = beta_cruise/(n_pr*alpha_cruise*P0_ICE_P0)*(2/(sigma*rho_SL)*beta_cruise*W0_S)^(1/2)*2*(Cd0_cruise*K^3)^(1/4)*(u^4+1)/(2*u)

% % Plotting out P0_W0_cruise as a function of height and W0_S
% P0_W0_cruise_array = []
% h_array = linspace(200,5000,25)
% W0_S_array = linspace(1000, 8000, 25)
% for iH = 1:length(h_array)
%     for iW0_S = 1:length(W0_S_array)
%         W3_W2_temp = double(vpasolve(subs(eq0,[h,W0_S],[h_array(iH),W0_S_array(iW0_S)]),W3_W2));
%         P0_W0_cruise_array(iW0_S,iH) = double(subs(P0_W0_cruise,[h,W0_S,W3_W2],[h_array(iH),W0_S_array(iW0_S),W3_W2_temp]));
%     end
% end
% figure
% [h_array, W0_S_array] = meshgrid(h_array, W0_S_array)
% surf(h_array,W0_S_array,P0_W0_cruise_array)
% xlabel('h (m)')
% ylabel('$\frac{W_0}{S}$')
% zlabel('$\frac{P_0}{W_0}$')

% Cruise altitude
h_cruise = 1000*0.3048
P0_W0_cruise_opt = subs(P0_W0_cruise,h,h_cruise)

W0_S_cruise = linspace(lower_W0_S,upper_W0_S,resolution);
W3_W2_cruise = []
P0_W0_cruise_opt_array = []
Cl_cruise_array = []
for iW0_S = 1:length(W0_S_cruise)
    W3_W2_cruise(iW0_S) = double(vpasolve(subs(eq0,[h,W0_S],[h_cruise,W0_S_cruise(iW0_S)]),W3_W2));
    P0_W0_cruise_opt_array(iW0_S)= double(subs(P0_W0_cruise_opt,[W0_S,W3_W2],[W0_S_cruise(iW0_S), W3_W2_cruise(iW0_S)]));
    Cl_cruise_array(iW0_S) = double(subs(Cl_cruise,[h,W0_S,W3_W2],[h_cruise,W0_S_cruise(iW0_S),W3_W2_cruise(iW0_S)]));
%     Cl_cruise = double(subs(Cl_cruise, W0_S, W0_S_cruise))
end


plot(W0_S_cruise,P0_W0_cruise_opt_array)

% Service ceiling
% Given a service ceiling of 25000
% Assuming a Vvmax of 100ft/min
service_ceiling = 25000 %ft
service_ceiling = service_ceiling*0.3048
sigma_SC = double(subs(sigma,h,service_ceiling))
alpha_climb = 0.9*sigma_SC^0.8 % Set 6 Slide 5
beta_climb = W1_W0*W2_W1

Vv_max_req = 100 % ft/min
Vv_max_req = Vv_max_req*0.3048/60 %m/s

P0_W0_climb = beta_climb/(alpha_climb)*Vv_max_req/n_pr+2*beta_climb^(3/2)/(n_pr*alpha_climb*sigma_SC*rho_SL)*(K/(3*Cd0_climb))^(1/2)*W0_S^(1/2)*1.155*beta_climb^(1/2)/L_D_star
W0_S_climb = linspace(lower_W0_S,upper_W0_S,resolution);
P0_W0_climb = double(subs(P0_W0_climb,W0_S,W0_S_climb));

plot(W0_S_climb,P0_W0_climb)

% OEI Climb Gradient
k_to = 1.2
h_OEI_climb = 0
sigma_OEI_climb = double(subs(sigma,h,h_OEI_climb))
beta_OEI_climb  = 1
alpha_OEI_climb = sigma_OEI_climb^0.8
sin_theta_req = 0.024
Cd0_OEI_climb = Cd0
eq_1 = sin_theta_req == alpha_OEI_climb*n_pr/beta_OEI_climb*P0_W0*1/k_to*(2*beta_OEI_climb/(sigma_OEI_climb*rho_SL)*W0_S*1/Cl_max_to)^(-1/2)-1/beta_OEI_climb*1/2*sigma_OEI_climb*rho_SL*(2*beta_OEI_climb/(sigma_OEI_climb*rho_SL)*W0_S*1/Cl_max_to)*k_to^2*W0_S^-1*(Cd0_OEI_climb+K*Cl_max_to/k_to^2)
P0_W0_OEI_climb = solve(eq_1,P0_W0)
W0_S_OEI_climb = linspace(lower_W0_S,upper_W0_S,resolution);

P0_W0_OEI_climb = double(subs(P0_W0_OEI_climb,W0_S,W0_S_OEI_climb));

plot(W0_S_OEI_climb,P0_W0_OEI_climb)


xlim([lower_W0_S-100,upper_W0_S + 100])
ylim([0,90])
% title(num2str([P0_ICE_P0,P0_elec_P0,P0_elec_prop_to,P0_elec_prop_climb],'Constraint plot for $\frac{P_{0_ICE}}{P_0}$=%.2f,$\frac{P_{0_e}}{P_0}$=%.2f,$P_{e_to}$=%.2f,$P_{e_climb}$=%.2f'))
title(num2str([P0_ICE_P0,P0_elec_P0],'Constraint plot for $P_{0_{ICE}}/P_0$=%.2f,$P_{0_e}/P_0$=%.2f'))

% title(num2str([1],'$W_0$ = %.1i'))
xlabel('$\frac{W_0}{S}$')
ylabel('$\frac{P_0}{W_0}$')
grid on
legend('Landing approach speed','TOFL','Cruise','Service Ceiling','OEI Climb Gradient')

%% Calculating preliminary weight using an expected P0/W0 and chosen W0/S
% Used for initial constraint analysis
% First run
% W0_S_chosen = 3050.51 %Pa
% For 6300kg payload (from prelim report)
W0_S_chosen = 2995

% First run
% P0_W0_aim = 31.9646 %W/N
% For 6300kg payload (from prelim report)
P0_W0_aim = 27.3599

% Cl at cruise (may be different from Cl_star)
sigma_cruise = double(subs(sigma,h,h_cruise))
W3_W2_chosen = double(vpasolve(subs(eq0,[h,W0_S],[h_cruise,W0_S_chosen])))
beta_cruise = double(subs(beta_cruise,[h,W3_W2,W0_S],[h_cruise,W3_W2_chosen,W0_S_chosen]))
Cl_cruise_chosen = 2/(sigma_cruise*rho_SL*V_cruise_req^2)*beta_cruise*W0_S_chosen

% Lift on drag at cruise
L_D_cruise_chosen = 1/(Cd0_cruise/Cl_cruise_chosen+K*Cl_cruise_chosen)

% Weight fraction from takeoff to landing
W4_W0_chosen = W1_W0*W2_W1*W3_W2_chosen*W4_W3

% Mass of fuel
fuel_reserve = 0.05
Wf_W0_chosen_harmonic = (1+fuel_reserve)*(1-W4_W0_chosen)

% Getting a more accurate estimate of the weight
M_payload = 6300 %kg
W_payload = M_payload*g

% Battery power output
battery_PO = 1.33 %kW/kg
% Motor efficiency
n_e = 0.95
W_batteries_W0_aim = P0_elec_P0*1/(n_e*battery_PO*1000)*P0_W0_aim

iterate = 10
W0_aim = 20000
for i = 1:iterate
    M0_aim = W0_aim/g
    We_W0_aim = 0.92*M0_aim^-0.05
    %Update weight
    W0_aim = W_payload/(1-Wf_W0_chosen_harmonic-We_W0_aim-W_batteries_W0_aim);
end
M0_aim = W0_aim/g

%% Read off the constraint curve and calculate performance parameters
% ICE Engine

% For Rolls Royce AE2100A (for hybrid)
% Power of individual ICE
P0_ICE_pw_chosen = 3096 %kW
% RPM of ICE
RPM_ICE_chosen = 1100 % RPM
% Mass of ICE
M_ICE_chosen = 715.77 %kg
% Number of ICE per wing
no_ICE_pw = 1

% % For Rolls Royce AE2100D2 (for fuel only)
% % Power of individual ICE
% P0_ICE_pw_chosen = 3458 %kW
% % RPM of ICE
% RPM_ICE_chosen = 1020.7 % RPM
% % Mass of ICE
% M_ICE_chosen = 805.49 %kg
% % Number of ICE per wing
% no_ICE_pw = 1

% P0_ICE_chosen = 2*P0_ICE_chosen*740.5 % W

P0_ICE_chosen = no_ICE_pw*2*P0_ICE_pw_chosen*1000 % W

% Electric motor

% For Pipistrel (for hybrid)
% Power of individual motor
P0_elec_pw_chosen = 57.6% kW % Pipistrel 
% RPM of motor
RPM_elec_chosen = 2500 % RPM
% Mass of motor
M_motor_chosen = 22.7 %kg
% Number of motor per wing
no_motor_pw = 6

% % No electric motor
% % Power of individual motor
% P0_elec_pw_chosen = 0% kW % Pipistrel 
% % RPM of motor
% RPM_elec_chosen = 2500 % RPM
% % Mass of motor
% M_motor_chosen = 0 %kg
% % Number of motor per wing
% no_motor_pw = 6

P0_elec_chosen = no_motor_pw*2*P0_elec_pw_chosen*1000 % W

P0_chosen = P0_ICE_chosen + P0_elec_chosen

P0_elec_P0_chosen = P0_elec_chosen/P0_chosen
P0_ICE_P0_chosen = P0_ICE_chosen/P0_chosen

% RPM chosen (take the lower of the two)
RPM_chosen = min([RPM_ICE_chosen,RPM_elec_chosen]) % RPM

% For a given percentage usage of the electric motor, aim for the below
% parameters
% P0_W0_aim = 34.9646 %W/N
% P0_aim = P0_W0_aim*W0
% P0_ICE_aim = P0_ICE_P0*P0_aim
% P0_ICE_aim_per_wing = P0_ICE_aim/2
% P0_elec_aim = P0_elec_P0*P0_aim
% P0_elec_aim_per_wing = P0_elec_aim/2

% Battery mass

m_batteries_chosen = P0_elec_chosen/(n_e)/(battery_PO*1000)
W_batteries_chosen = m_batteries_chosen*g

% Weight of plane at MTOW
iterate = 10
W0_chosen = W0_aim
for i = 1:iterate
    M0_chosen = W0_chosen/g
    W_batteries_W0_chosen = W_batteries_chosen/W0_chosen
    We_W0_chosen = 0.92*M0_chosen^-0.05
    %Update weight
    W0_chosen = W_payload/(1-Wf_W0_chosen_harmonic-We_W0_chosen-W_batteries_W0_chosen);
    % Update parameters
    S_chosen = W0_S_chosen^-1*W0_chosen
    P0_W0_chosen = P0_chosen/W0_chosen
end
M0_chosen = W0_chosen/g

% Fuel mass
Wf_chosen_harmonic = Wf_W0_chosen_harmonic*W0_chosen
Mf_chosen_harmonic = Wf_chosen_harmonic/g

% From Qunhy's VSP
Mf_chosen_max = 3500
Wf_chosen_max = Mf_chosen_max*g

% Empty mass
Me_chosen = 0.92*M0_chosen^-0.05*W0_chosen/g
We_chosen = Me_chosen*g

% Span
b_chosen = (A*S_chosen)^(1/2)

% Propeller sizing (need to ensure V_helical < 240)
% For four propellers
n_propellers = 4;
if n_propellers == 4
    D_p = 0.49*(P0_chosen/1000/2)^0.25
end

V_helical = ((pi*RPM_chosen/60*D_p)^2+V_cruise_req^2)^(1/2)

D_p_chosen = 3.1 %m
V_helical_chosen = ((pi*RPM_chosen/60*D_p_chosen)^2+V_cruise_req^2)^(1/2)

% Fuselage length
l_fuselage_chosen = 0.169*M0_chosen^0.51 %m (Correlations for twin turboprops by Raymer)

% Tail moment arms

% Horizontal tail
L_HT = 0.525*l_fuselage_chosen %m (Rule of thumb by Raymer)

% Vertical tail
L_VT = L_HT % Assuming vertical tail length is the same as horizontal

% Tail surface areas

% Geometric mean chord
c_gm = S_chosen/b_chosen

% Wing MAC 
c_mac_w = c_gm % Assume it is close to geometric chord

% Correlations for horizontal and vertical tail by Raymer
c_HT = 0.9 
c_VT = 0.08

% Horizontal tail surface area
S_HT = c_HT*c_mac_w*S_chosen/L_HT %m^2

% Vertical tail surface area
S_VT = c_VT*b_chosen*S_chosen/L_VT %m^2

% Allow for a 5% reduction so that vertical taile puts horizontal tail out
% of up-wash
S_HT = 0.95*S_HT
S_VT = 0.95*S_VT



% Harmonic Range (Take off at MPL and MTOW)
recharge_reserve_harmonic = 0.05
reserve_fuel_harmonic = Wf_chosen_harmonic*fuel_reserve
recharge_fuel_harmonic = Wf_chosen_harmonic*recharge_reserve_harmonic
Wf_available_harmonic = Wf_chosen_harmonic-reserve_fuel_harmonic-recharge_fuel_harmonic
W0_harmonic = Wf_chosen_harmonic + W_payload + We_chosen + W_batteries_chosen
W1_harmonic = W1_W0*W0_harmonic
Wf_TO_harmonic = W0_harmonic-W1_harmonic
W2_harmonic = W2_W1*W1_harmonic
Wf_climb_harmonic = W1_harmonic - W2_harmonic
% Notes to determine W3_harmonic
% fuel_cruise_harmonic = avaialable_fuel_harmonic - fuel_TO_harmonic -
% fuel_climb_harmonic - fuel_land_harmonic
% fuel_land_harmonic = W4_harmonic - W3_harmonic
% fuel_cruise_harmonic = avaialable_fuel_harmonic - fuel_TO_harmonic -
% fuel_climb_harmonic - (W4_harmonic - W3_harmonic)
% W3_harmonic = W2_harmonic - (avaialable_fuel_harmonic - fuel_TO_harmonic -
% fuel_climb_harmonic - (W4_harmonic - W3_harmonic))
% W4_harmonic = W4_W3*W3_harmonic
% After manipulating...
W3_harmonic = (W2_harmonic-Wf_available_harmonic+Wf_TO_harmonic + Wf_climb_harmonic)/(2-W4_W3)
R_h_chosen = n_pr/(g*cp_cruise)*L_D_cruise_chosen*log(W2_harmonic/W3_harmonic) %km
R_h_chosen = R_h_chosen/1000/1.852

% Ferry Range (Max fuel capacity and no payload)
W0_ferry = We_chosen + Wf_chosen_max + W_batteries_chosen
Wf_total = Wf_chosen_max
recharge_reserve_ferry = 0.05
reserve_fuel_ferry = Wf_total*0.05
recharge_fuel_ferry = Wf_total*recharge_reserve_ferry
Wf_available_ferry = Wf_total - reserve_fuel_ferry-recharge_fuel_ferry
W1_ferry = W1_W0*W0_ferry
Wf_TO_ferry = W0_ferry-W1_ferry
W2_ferry = W2_W1*W1_ferry
Wf_climb_ferry = W1_ferry - W2_ferry
W3_ferry = (W2_ferry-Wf_available_ferry+Wf_TO_ferry + Wf_climb_ferry)/(2-W4_W3)
R_f_chosen = n_pr/(g*cp_cruise)*L_D_cruise_chosen*log(W2_ferry/W3_ferry) %km
R_f_chosen = R_f_chosen/1000/1.852

% Max Range (Max fuel capacity and MTOW)
W0_max = W0_chosen
Wf_total = Wf_chosen_max
recharge_reserve_max = 0.05
reserve_fuel_max = Wf_total*0.05
recharge_fuel_max = Wf_total*recharge_reserve_max
Wf_available_max = Wf_total - reserve_fuel_max-recharge_fuel_max
W1_max = W1_W0*W0_max
Wf_TO_max = W0_max-W1_max
W2_max = W2_W1*W1_max
Wf_climb_max = W1_max - W2_max
W3_max = (W2_max-Wf_available_max+Wf_TO_max + Wf_climb_max)/(2-W4_W3)
R_max_chosen = n_pr/(g*cp_cruise)*L_D_cruise_chosen*log(W2_max/W3_max) %km
R_max_chosen = R_max_chosen/1000/1.852



% TOFL
TOP_chosen = P0_W0_chosen^-1*W0_S_chosen*1/sigma_to*1/Cl_max_to
TOFL_chosen = 11.8*TOP_chosen + 0.255*TOP_chosen^2

% Service ceiling
sigma_climb_chosen = sigma
alpha_climb_chosen = 0.9*sigma_climb_chosen^0.8
beta_climb_chosen = W1_W0*W2_W1

eq_2 = 100*0.3048/60 == n_pr*alpha_climb_chosen/beta_climb_chosen*P0_W0_chosen-2/(sigma_climb_chosen*rho_SL)*(K/(3*Cd0_climb))^(1/2)*(W0_S_chosen)^(1/2)*(1.155*beta_climb_chosen^(1/2)/L_D_star)
assume(h,'real')
service_ceiling_chosen = double(vpasolve(eq_2,h))
service_ceiling_chosen = service_ceiling_chosen/0.3048

% Absolute ceiling
eq_3 = 0 == n_pr*alpha_climb_chosen/beta_climb_chosen*P0_W0_chosen-2/(sigma_climb_chosen*rho_SL)*(K/(3*Cd0_climb))^(1/2)*(W0_S_chosen)^(1/2)*(1.155*beta_climb_chosen^(1/2)/L_D_star)
assume(h,'real')
absolute_ceiling = double(vpasolve(eq_3,h))
absolute_ceiling = absolute_ceiling/0.3048

% Summarizing results
disp(num2str(P0_elec_P0_chosen*100,'Specification for %.2f percent electrical contribution:'))
disp(num2str([P0_elec_pw_chosen,no_motor_pw],'TO Power per motor = %.2f kW (x %.1i per wing)'))
try
    disp(num2str(M_motor_chosen, 'Motor mass = %.2f kg'))
catch
end
disp(num2str(RPM_elec_chosen, 'Motor RPM = %.2f'))
disp(num2str([P0_ICE_chosen,no_ICE_pw],'TO Power per ICE = %.2f kW (x %.1i per wing)'))
try
    disp(num2str(M_ICE_chosen, 'ICE mass = %.2f kg'))
catch
end

disp(num2str(RPM_ICE_chosen, 'ICE RPM = %.2f'))
disp(num2str(P0_elec_chosen/1000,'Total electrical power = %.2f kW'))
disp(num2str(P0_ICE_chosen/1000,'Total ICE power = %.2f kW'))
disp(num2str(M_motor_chosen*no_motor_pw*2,'Total motor mass = %.2f kg'))
disp(num2str(M_ICE_chosen*no_ICE_pw*2,'Total ICE mass = %.2f kg'))
disp(num2str(M0_chosen,'Maximum takeoff mass = %.2f kg'))
disp(num2str(Me_chosen,'Empty mass = %.2f kg'))
disp(num2str(Mf_chosen,'Fuel mass = %.2f kg (for harmonic range)'))
disp(num2str(m_batteries_chosen,'Battery mass = %.2f kg'))
disp(num2str(l_fuselage_chosen,'Fuselage length = %.2f m'))
disp(num2str(A,'Aspect Ratio = %.2f'))
disp(num2str(S_chosen,'Wing area = %.2f m^2'))
disp(num2str(b_chosen,'Wing span = %.2f m'))
disp(num2str(L_HT,'Horizontal tail moment arm = %.2f m'))
disp(num2str(S_HT,'Horizontal tail area = %.2f m^2'))
disp(num2str(L_VT,'Vertical tail moment arm = %.2f m'))
disp(num2str(S_VT,'Vertical tail area = %.2f m^2'))
disp(num2str([D_p_chosen, n_propellers],'Propeller diameter = %.2f m (x %.1i blades)'))
disp(num2str(service_ceiling_chosen,'Service ceiling = %.1f ft'))
disp(num2str(absolute_ceiling,'Absolute ceiling = %.1f ft'))
disp(num2str(R_h_chosen,'Harmonic range = %.1f nm'))
disp(num2str(R_max_chosen,'Maximum range = %.1f nm'))
disp(num2str(R_f_chosen,'Ferry range = %1f nm'))
disp(num2str(TOFL_chosen,'TOFL = %1f m'))




%% Non-constraint analysis (because Jayson says that constraint analysis is overkill but can be used for the appendices)
% For takeoff

% Using preliminary estimated mass
M0 = 1.949591926490346e+04 %kg
%Gravity
g = 9.8065
% Weight
W0 = M0*g
% Using estimates of W0_S from Raymer
W0_S = 3000
S = W0_S^-1*W0
% Sigma at SL
sigma_to = 1
% Cl max at TO
Cl_max_to = 2.0
syms P0
% Takeoff parameter
TOP = W0/P0*W0_S*1/sigma_to*1/Cl_max_to
%TOFL for two engine
TOFL = 11.8*TOP + 0.255*TOP^2
%TOFL to aim at
% For STOL
TOFL_r = 1000 
% Solving for  P0
assume(P0,'positive')
P0_to = double(solve(TOFL == TOFL_r,P0))
% Power distribution
% Motor
P0_elec_P0_to = 0.2728
% ICE
P0_ICE_P0_to = 1-P0_elec_P0_to
% Power motor for to
P0_elec_to = P0_elec_P0_to*P0_to
% Power ICE for to
P0_ICE_to = P0_ICE_P0_to*P0_to

%Battery power output
batt_max_P_out = 1.33 %kW/kg
%Motor efficiency
n_e = 0.95
%Mass of batteries
m_batt=  P0_elec_to/n_e/(batt_max_P_out*1000) % This is the mass of batteries required to drive the motor at full power
% Battery capacity
e_batt = 0.333 %kWh/kg
E_batt = m_batt*e_batt %kWh
E_batt = E_batt*1000*3600 %J


% Takeoff factor on stall speed
k_to = 1.2
% Lift at takeoff
Cl_2 = Cl_max_to/k_to^2
% Sea level density
rho_SL = 1.225
% Velocity at takeoff
V_2 = (2/rho_SL*W0_S*1/Cl_2)^(1/2) %m/s
% Propeller blades
K_p = 0.49
% Blade diameter and area
D_p = K_p*(P0_to/1000)^(1/4)
A_p = (pi*D_p^2/4)
% Vary velocity
syms V
% First segment of thrust
T_seg_1 = (2*P0_to^2*rho_SL*A_p)^(1/3)
% Propulsive efficiency
n_p = 0.85
% Second segment of thrust
T_seg_2 = (n_p*n_e*P0_elec_to+n_p*P0_ICE_to)/V
% Intersecting velocity between T_seg_1 and 2
V_int = double(solve(T_seg_1==T_seg_2,V))
% Integrate the thrust over V = 0 to V = V2
T_integrated = double(int(T_seg_1,V,[0,V_int]) + int(T_seg_2,V,[V_int,V_2]))
% Average acceleration
a_average = T_integrated/M0*1/V_2
% Reduce
r_t = 0.8
a_average = 0.8*a_average
% Roll runway
s_run = V_2^2/(2*a_average)
% Clearance height
h_to = 35 %ft
h_to = h_to*0.3048
T_2 = double(subs(T_seg_2,V,V_2))
% Aerodynamic parameters
Cd0_to = 0.02 % Torenbeek
e = 0.7
% Aspect Ratio
A = 12
K = 1/(pi*A*e)
% Drag at end of takeoff
Cd_2 = Cd0_to+K*Cl_2^2
% Weight factor
beta_to = 0.97
% Takeoff climb angle
tan_gamma_2 = (T_2/(beta_to*W0)-Cd_2/Cl_2)
gamma_2 = atan(tan_gamma_2)*180/pi
s_air = h_to/tan_gamma_2
% Takeoff
s_to = s_run + s_air
% Adding 15% margin
s_to = 1.15*s_to

% Time taken
t_run = V_2/a_average
t_air = h_to/(V_2*sin(gamma_2*pi/180))
t_to = t_run + t_air
% Voltage of motor
V_motor = 750 %V
% Battery capacity used
E_to = t_to*P0_elec_to
% Battery capacity left
E_left = E_batt-E_to


% Cruise

% Alitude density effect
syms h
T = -6.5*10^-3*h + 288.16 %K
sigma = (T/288.16)^(-g/(-6.5*10^-3*287)-1)
alpha = 0.7*sigma^0.8

% Required harmonic range
R_h = 700 %nm
R_h = R_h*1.852*1000 %m
% Cruise speed required
V_cr_req = 300 %KTAS
V_cr_req = V_cr_req/1.94384
% Zero drag
Cd0_cr = 0.02
% Max lift to drag ratio
L_D_star = 1/(4*Cd0_cr*K)^(1/2)
% Lift to drag ratio
Cl_cr = 2/(sigma*rho_SL*V_cr_req^2)*W0_S
Cd_cr = Cd0_cr + K*Cl_cr^2
L_D_cr = Cl_cr/Cd_cr
% PSFC
cp_cr = 0.085*10^-6 % kg/Ws
cp_cr = 0.7*cp_cr % Reduce due to parallel
% Gravity
g = 9.8065


W3_W2 = exp(-R_h*g*cp_cr/(n_p*L_D_cr))
% Cruise weight fraction
W1_W0 = 0.97
W2_W1 = 0.985
beta_cr = W1_W0*W2_W1*(W3_W2+1)/2

% Finding the optimal cruise altitude attempt 2

% Max lift to drag speed
Cl_star = (Cd0_cr/K)^(1/2)
V_star = (2/(sigma*rho_SL)*W0_S)^(1/2)*(1/Cl_star)^(1/2)

% Percent electricity and ICE
P0_elec_P0 = 0.1
P0_ICE_P0 = 1-P0_elec_P0

u = V_cr_req/V_star
P0_cr = beta_cr/(n_p*alpha*P0_ICE_P0)*(2/(sigma*rho_SL)*W0_S)^(1/2)*2*(Cd0_cr*K^3)^(1/4)*(u^4+1)/(2*u)*W0
figure
fplot(P0_cr)
xlim([0 25000*0.3048])


% Try to get min power and optimal altitude
% Derive P0_cr
dP0_cr_dV = diff(P0_cr,h)
% Optimal height
assume(h,'positive')
h_cr = double(solve(dP0_cr_dV == 0,h))
% Optimal power
P0_cr_opt = double(subs(P0_cr,h,h_cr))
% Alpha
alpha_cr_opt = double(subs(alpha,h,h_cr))


if length(P0_cr_opt) == 0 || h_cr < 1000*0.3048
    % Cruise altitude
    h_cr = 1000 % ft
    h_cr = h_cr*0.3048
    P0_cr_opt = double(subs(P0_cr,h,h_cr))
end
% Service ceiling
h_sc = 25000 %ft
h_sc = h_sc*0.3048
P0_sc = double(subs(P0_cr,h,h_sc))

% % If cruise is run on fuel then calculate the percentage of electricity
% % used for takeoff
% P0_ICE_P0_opt = P0_cr_opt/P0_to
% P0_elec_P0_opt = 1-P0_ICE_P0_opt

% % Range 
sigma_cr_opt = double(subs(sigma,h,h_cr))
W3_W2_cr_opt = double(subs(W3_W2,h,h_cr))
beta_cr_opt = subs(beta_cr,h,h_cr)
Cl_cr_opt = 2/(rho_SL*sigma_cr_opt*V_cr_req^2)*beta_cr_opt*W0_S
Cd_cr_opt = Cd0_cr + K*Cl_cr_opt^2
L_D_cr_opt = Cl_opt/Cd_opt
R_cr = double(n_p/(g*cp_cr)*L_D_cr_opt*log(W3_W2_cr_opt^-1)) %m
R_cr_nm = R_cr/1.852/1000






