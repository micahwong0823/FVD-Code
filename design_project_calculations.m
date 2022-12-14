%% Initial weight estimates (skip this for now, start at constraint analysis)
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

% W0_S axis when plotting
lower_W0_S = 1000
upper_W0_S = 8000
resolution = 100

figure_position = [100 100 1000 800]
constraint_plot = figure
constraint_plot.Position = figure_position

% Shared parameters
A = 11.5

g = 9.8065 %m/s^2

rho_SL = 1.225
e = 0.7 % Set 6 Slide 2

K = 1/(pi*A*e)
n_pr = 0.85
CD0 = 0.02 % Torenbeek

% Reduction to cp
HEP_cp_reduction = 0.3 % If hybrid parallel
% HEP_cp_reduction = 0 % If fuel only

k_land = 1.3
k_to = 1.2


service_ceiling = 25000 % ft
service_ceiling = service_ceiling*0.3048

fuel_reserve = 0.05

% Battery power output
battery_PO = 1.33 %kW/kg

cp_cruise = 0.085*10^-6 % kg/Ws
cp_cruise = cp_cruise*(1-HEP_cp_reduction)
cp_loiter = 0.101*10^-6 %kg/Ws
cp_loiter = cp_loiter*(1-HEP_cp_reduction)

cp_cruise_constraint = cp_cruise
cp_loiter_constraint = cp_loiter

% Analytical curve for sigma
syms h
T = -6.5*10^-3*h + 288.16 %K
sigma = (T/288.16)^(-g/(-6.5*10^-3*287)-1)

%Hybridization
P0_elec_P0_constraint = 0.1 
P0_ICE_P0_constraint = 1-P0_elec_P0_constraint

% Initial parameters
W1_W0_constraint = 0.97 % Note that if we are using electricity for this section then W1_W0 = 1
W2_W1_constraint = 0.985
W4_W3_constraint = 0.995

% Apply electric proportion
W1_W0_constraint = 1-P0_ICE_P0_constraint*(1-W1_W0_constraint)
W1_W0_constraint = 1-P0_ICE_P0_constraint*(1-W1_W0_constraint)
W1_W0_constraint = 1-P0_ICE_P0_constraint*(1-W1_W0_constraint)

% Note We_W0 = 0.92*M0^-0.05


CD0_to_constraint = CD0
CD0_climb_constraint = CD0
CD0_cruise_constraint = CD0
CD0_land_constraint = CD0

% A = 12 % Qunhy's assumption


e_to_constraint = e
e_climb_constraint = e
e_cruise_constraint = e
e_landing_constraint = e


K_to_constraint = K
K_climb_constraint = K
K_cruise_constraint = K
K_landing_constraint = K

CL_max_clean_constraint = 1.65
CL_max_climb_constraint = 1.65
CL_max_to_constraint = 1.95
CL_max_land_constraint = 2.3 % Torenbeek ?

n_pr_cruise_constraint = n_pr
n_pr_climb_constraint = n_pr
n_pr_to_constraint = n_pr
n_pr_landing_constraint = n_pr


% Landing approach speed 
% Assume approaching at Vapp = 1.3*V_stall 
% Also assume landing at SL
% Assume beta_land = 1
k_land_constraint = 1.3
V_app_req = 77 %m/s https://skybrary.aero/articles/approach-speed-categorisation
sigma_land_constraint = 1
beta_land_constraint = 1 % Not sure about this assumption
rho_land_constraint = rho_SL
% From set 6 slides (Hugh)
W0_S_land_constraint = (V_app_req/k_land_constraint)^2*(sigma_land_constraint*rho_land_constraint/2)*1/beta_land_constraint*CL_max_land_constraint % Pa

xline(W0_S_land_constraint)
hold on

% Take off field length
% Assume taking off at SL
% TOFL = 1000m (for short take off)
syms W0_S P0_W0
assume(W0_S,'positive')
assume(P0_W0,'positive')

sigma_to_constraint = 1

% From set 6 slides (Hugh)
TOP_constraint = W0_S*(P0_W0)^-1*1/sigma_to_constraint*1/CL_max_to_constraint
TOFL_constraint = 11.8*TOP_constraint+0.255*TOP_constraint^2 % For two engine prop
TOFL_req = 1000

W0_S_to_constraint = linspace(lower_W0_S,upper_W0_S,resolution);
P0_W0_to_constraint = []
for iW0_S = 1:length(W0_S_to_constraint)
    P0_W0_to_temp = solve(TOFL_req == subs(TOFL_constraint,W0_S,W0_S_to_constraint(iW0_S)),P0_W0);
    P0_W0_to_constraint(iW0_S) = double(P0_W0_to_temp);
end

plot(W0_S_to_constraint,P0_W0_to_constraint)

% Cruise speed

% Alitude density effect
syms W3_W2
alpha_cruise_constraint = 0.7*sigma^0.8

% Required harmonic range
R_h = 700 %nm
R_h = R_h*1.852*1000 %m

% Cruise speed required
V_cruise_req = 300 %KTAS
V_cruise_req = V_cruise_req/1.94384

% Max lift to drag ratio
L_D_star_cruise_constraint = 1/(4*CD0_cruise_constraint*K_cruise_constraint)^(1/2)

% Cruise weight fraction
beta_cruise_constraint = W1_W0_constraint*W2_W1_constraint*(W3_W2+1)/2

% Lift to drag ratio
CL_cruise_constraint = 2/(sigma*rho_SL*V_cruise_req^2)*beta_cruise_constraint*W0_S
CD_cruise_constraint = CD0_cruise_constraint + K*CL_cruise_constraint^2
L_D_cruise_constraint = CL_cruise_constraint/CD_cruise_constraint

% Solving cruise performance equation (set 13 slides)
eq0 = W3_W2 == exp(-R_h*g*cp_cruise_constraint/(n_pr_cruise_constraint*L_D_cruise_constraint))
% Max lift to drag speed
CL_star_cruise_constraint = (CD0_cruise_constraint/K_cruise_constraint)^(1/2)
V_star_cruise_constraint = (2/(sigma*rho_SL)*beta_cruise_constraint*W0_S)^(1/2)*(1/CL_star_cruise_constraint)^(1/2)
u = V_cruise_req/V_star_cruise_constraint
P0_W0_cruise_constraint = beta_cruise_constraint/(n_pr_cruise_constraint*alpha_cruise_constraint*P0_ICE_P0_constraint)*(2/(sigma*rho_SL)*beta_cruise_constraint*W0_S)^(1/2)*2*(CD0_cruise_constraint*K_cruise_constraint^3)^(1/4)*(u^4+1)/(2*u)

% Cruise altitude
h_cruise_constraint = 22000*0.3048
P0_W0_cruise_opt_constraint = subs(P0_W0_cruise_constraint,h,h_cruise_constraint)

W0_S_cruise_constraint = linspace(lower_W0_S,upper_W0_S,resolution);
W3_W2_cruise_constraint = []
P0_W0_cruise_opt_array_constraint = []
CL_cruise_array_constraint = []
% Solving for W3_W2, P0_W0 and CL for a bunch of W0_S
for iW0_S = 1:length(W0_S_cruise_constraint)
    W3_W2_cruise_constraint(iW0_S) = double(vpasolve(subs(eq0,[h,W0_S],[h_cruise_constraint,W0_S_cruise_constraint(iW0_S)]),W3_W2));
    P0_W0_cruise_opt_array_constraint(iW0_S)= double(subs(P0_W0_cruise_opt_constraint,[W0_S,W3_W2],[W0_S_cruise_constraint(iW0_S), W3_W2_cruise_constraint(iW0_S)]));
    CL_cruise_array_constraint(iW0_S) = double(subs(CL_cruise_constraint,[h,W0_S,W3_W2],[h_cruise_constraint,W0_S_cruise_constraint(iW0_S),W3_W2_cruise_constraint(iW0_S)]));
end


plot(W0_S_cruise_constraint,P0_W0_cruise_opt_array_constraint)


% Service ceiling
% Given a service ceiling of 25000
% Assuming a Vvmax of 100ft/min
service_ceiling_constraint = service_ceiling %ft
sigma_SC_constraint = double(subs(sigma,h,service_ceiling_constraint))
alpha_climb_constraint = 0.9*sigma_SC_constraint^0.8 % Set 6 Slide 5
beta_climb_constraint = W1_W0_constraint*W2_W1_constraint

Vv_max_req = 100 % ft/min
Vv_max_req = Vv_max_req*0.3048/60 %m/s
L_D_star_SC_constraint = L_D_star_cruise_constraint
% Set 6 slides (Hugh)
P0_W0_climb_constraint = beta_climb_constraint/(alpha_climb_constraint)*Vv_max_req/n_pr_climb_constraint+2*beta_climb_constraint^(3/2)/(n_pr_climb_constraint*alpha_climb_constraint*sigma_SC_constraint*rho_SL)*(K_climb_constraint/(3*CD0_climb_constraint))^(1/2)*W0_S^(1/2)*1.155*beta_climb_constraint^(1/2)/L_D_star_SC_constraint
W0_S_climb_constraint = linspace(lower_W0_S,upper_W0_S,resolution);
P0_W0_climb_constraint = double(subs(P0_W0_climb_constraint,W0_S,W0_S_climb_constraint));

plot(W0_S_climb_constraint,P0_W0_climb_constraint)

% OEI Climb Gradient
k_to_constraint = 1.2
h_OEI_climb_constraint = 0
sigma_OEI_climb_constraint = double(subs(sigma,h,h_OEI_climb_constraint))
beta_OEI_climb_constraint  = 1
alpha_OEI_climb_constraint = sigma_OEI_climb_constraint^0.8
sin_theta_req = 0.024
CD0_OEI_climb_constraint = CD0
n_pr_OEI_climb_constraint = n_pr
K_OEI_climb_constraint = K
% Set 6 slides (Hugh)
eq_1 = sin_theta_req == alpha_OEI_climb_constraint*n_pr_OEI_climb_constraint/beta_OEI_climb_constraint*P0_W0*1/k_to_constraint*(2*beta_OEI_climb_constraint/(sigma_OEI_climb_constraint*rho_SL)*W0_S*1/CL_max_to_constraint)^(-1/2)-1/beta_OEI_climb_constraint*1/2*sigma_OEI_climb_constraint*rho_SL*(2*beta_OEI_climb_constraint/(sigma_OEI_climb_constraint*rho_SL)*W0_S*1/CL_max_to_constraint)*k_to_constraint^2*W0_S^-1*(CD0_OEI_climb_constraint+K_OEI_climb_constraint*CL_max_to_constraint/k_to_constraint^2)
P0_W0_OEI_climb_constraint = solve(eq_1,P0_W0)
W0_S_OEI_climb_constraint = linspace(lower_W0_S,upper_W0_S,resolution);

P0_W0_OEI_climb_constraint = double(subs(P0_W0_OEI_climb_constraint,W0_S,W0_S_OEI_climb_constraint));

plot(W0_S_OEI_climb_constraint,P0_W0_OEI_climb_constraint)


% Turning Constraint (assume this is done during cruise)
bank_angle_req = 30 % deg
bank_angle_req = bank_angle_req*pi/180
n_bank_req = cos(bank_angle_req)^-1
h_bank_constraint = h_cruise_constraint
beta_bank_constraint = W1_W0_constraint*W2_W1_constraint % Start of cruise
alpha_bank_constraint = 0.7*sigma^0.8
alpha_bank_constraint = double(subs(alpha_bank_constraint,h,h_bank_constraint))
sigma_bank_constraint = double(subs(sigma,h,h_bank_constraint))
q_bank_constraint = 1/2*sigma_bank_constraint*rho_SL*V_cruise_req^2
% Set 6 slides (Hugh)
P0_W0_bank_turn_constraint = V_cruise_req*beta_bank_constraint/(n_pr_cruise_constraint*alpha_bank_constraint)*(CD0_cruise_constraint*q_bank_constraint/(beta_bank_constraint*W0_S)+K_cruise_constraint*n_bank_req^2*beta_bank_constraint*W0_S/q_bank_constraint)

W0_S_bank_climb_constraint = linspace(lower_W0_S,upper_W0_S,resolution);
P0_W0_bank_turn_constraint = double(subs(P0_W0_bank_turn_constraint,W0_S,W0_S_OEI_climb_constraint));

plot(W0_S_bank_climb_constraint,P0_W0_bank_turn_constraint)


xlim([lower_W0_S-100,upper_W0_S + 100])
ylim([0,90])
title(num2str([P0_ICE_P0_constraint,P0_elec_P0_constraint],'Constraint plot for $P_{0_{ICE}}/P_0$=%.2f,$P_{0_e}/P_0$=%.2f'))

xlabel('$\frac{W_0}{S}$')
ylabel('$\frac{P_0}{W_0}$')
grid on
legend('Landing approach speed','TOFL','Cruise','Service Ceiling','OEI Climb Gradient','Bank Turn')

%% Calculating preliminary weight using an expected P0/W0 and chosen W0/S
% Following calculations are for values to aim at when choosing the actual
% aircraft parameters

% Used for initial constraint analysis
% For 6300kg payload (from prelim report)
W0_S_chosen = 3000

% For 6300kg payload (from prelim report)
P0_W0_aim = 27.3599

% CL at cruise (may be different from CL_star)
W3_W2_aim_constraint = double(vpasolve(subs(eq0,[h,W0_S],[h_cruise_constraint,W0_S_chosen]))) % Note that this is for harmonic range
beta_cruise_chosen_constraint = double(subs(beta_cruise_constraint,[h,W3_W2,W0_S],[h_cruise_constraint,W3_W2_aim_constraint,W0_S_chosen]))
sigma_cruise_constraint = double(subs(sigma,h,h_cruise_constraint))
CL_cruise_aim_constraint = 2/(sigma_cruise_constraint*rho_SL*V_cruise_req^2)*beta_cruise_chosen_constraint*W0_S_chosen
CL_cruise_aim_constraint = double(subs(CL_cruise_aim_constraint,[h,W0_S],[h_cruise_constraint,W0_S_chosen]))

% Lift on drag at cruise
L_D_cruise_aim_constraint = 1/(CD0_cruise_constraint/CL_cruise_aim_constraint+K_cruise_constraint*CL_cruise_aim_constraint)

% Weight fraction from takeoff to landing
W4_W0_aim_constraint = W1_W0_constraint*W2_W1_constraint*W3_W2_aim_constraint*W4_W3_constraint

% Mass of fuel
Wf_W0_aim_harmonic = (1+fuel_reserve)*(1-W4_W0_aim_constraint)

% Getting a more accurate estimate of the weight
M_payload_max = 6720 %kg
W_payload_max = M_payload_max*g


% Motor efficiency
n_e = 0.95
W_batteries_W0_aim = P0_elec_P0_constraint*1/(n_e*battery_PO*1000)*P0_W0_aim*g

iterate = 10
W0_aim = 20000
for i = 1:iterate
    M0_aim = W0_aim/g
    We_W0_aim = 0.92*M0_aim^-0.05
    %Update weight
    W0_aim = W_payload_max/(1-Wf_W0_aim_harmonic-We_W0_aim-W_batteries_W0_aim);
end
M0_aim = W0_aim/g
Wf_max_aim = Wf_W0_aim_harmonic*W0_aim
Mf_max_aim = Wf_max_aim/g
W_batteries_aim = W_batteries_W0_aim*W0_aim
M_batteries_aim = W_batteries_aim/g
We_aim = We_W0_aim*W0_aim
Me_aim = We_aim/g

%% Read off the constraint curve and calculate performance parameters
% ICE Engine

% For Rolls Royce AE2100A (for hybrid)
% Power of individual ICE
P0_ICE_pw_chosen = 3096 %kW
% P0_ICE_pw_chosen = 10000
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

% Total combustion power
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

% Total electrical power
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

% Battery mass (and other parameters)
M_batteries_chosen = 1419 %kg
W_batteries_chosen = M_batteries_chosen*g
battery_capacity = 0.333 % kWh/kg
battery_capacity = battery_capacity*1000*3600 %J/kg

% PO = power output
battery_max_PO = 1.33*10^3 % W/kg
battery_recharge_rate = 1.66*10^3 % W/kg

% Energy capacity
no_batteries = 4;
E_total_battery_chosen = M_batteries_chosen*battery_capacity
E_per_battery_chosen = E_total_battery_chosen/no_batteries;
E_battery_1_chosen = E_per_battery_chosen;
E_battery_2_chosen = E_per_battery_chosen;
E_battery_3_chosen = E_per_battery_chosen;
E_battery_4_chosen = E_per_battery_chosen;



% MTOW
% M0_MTOW_chosen = 23840
M0_MTOW_chosen = 22500
W0_MTOW_chosen = M0_MTOW_chosen*g

% OEW
Me_chosen = 12200
We_chosen = Me_chosen*g

% Fuel mass
Wf_aim_harmonic = Wf_W0_aim_harmonic*W0_MTOW_chosen
Mf_aim_harmonic= Wf_aim_harmonic/g

% From Qunhy's VSP
Mf_chosen_max = 6400 %kg
Wf_chosen_max = Mf_chosen_max*g

% Empty mass
Me_chosen = 0.92*M0_MTOW_chosen^-0.05*W0_MTOW_chosen/g
We_chosen = Me_chosen*g

% Span
A_aim = 11.5
S_chosen = 74
b_aim = (S_chosen*A_aim)^(1/2)
b_chosen = 28.531

% Propeller sizing (need to ensure V_helical < 240)
% For four propellers
n_propellers = 4;
if n_propellers == 4
    D_p_aim = 0.49*(P0_chosen/1000/2)^0.25
end

V_helical = ((pi*RPM_chosen/60*D_p_aim)^2+V_cruise_req^2)^(1/2)

D_p_chosen = 3.1 %m
V_helical_chosen = ((pi*RPM_chosen/60*D_p_chosen)^2+V_cruise_req^2)^(1/2)

% Fuselage length
l_fuselage_aim = 0.169*M0_MTOW_chosen^0.51 %m (Correlations for twin turboprops by Raymer)
l_fuselage_chosen = 27
% Tail moment arms

% Horizontal tail
L_HT_aim = 0.525*l_fuselage_aim %m (Rule of thumb by Raymer)

% Vertical tail
L_VT_aim = L_HT_aim % Assuming vertical tail length is the same as horizontal

% Tail surface areas

% Geometric mean chord
c_gm_aim = S_chosen/b_aim

% Wing MAC 
c_mac_w_aim = c_gm_aim % Assume it is close to geometric chord

% Correlations for horizontal and vertical tail by Raymer
c_HT = 0.9 
c_VT = 0.08

% Horizontal tail surface area
S_HT_aim = c_HT*c_mac_w_aim*S_chosen/L_HT_aim %m^2

% Vertical tail surface area
S_VT_aim = c_VT*b_aim*S_chosen/L_VT_aim %m^2

% Allow for a 5% reduction so that vertical taile puts horizontal tail out
% of up-wash
S_HT_aim = 0.95*S_HT_aim
S_VT_aim = 0.95*S_VT_aim

cp_to_chosen = cp_cruise
cp_climb_chosen = cp_cruise
cp_cruise_chosen = cp_cruise
cp_descent_chosen = cp_cruise
cp_landing_chosen = cp_cruise

K_to_chosen = K
K_climb_chosen = K
K_cruise_chosen = K
K_descent_chosen = K
K_landing_chosen = K

CD0_to_chosen = CD0
CD0_climb_chosen = CD0
CD0_cruise_chosen = CD0
CD0_descent_chosen = CD0
CD0_landing_chosen = CD0

n_pr_to_chosen = n_pr
n_pr_climb_chosen = n_pr
n_pr_cruise_chosen = n_pr
n_pr_descent_chosen = n_pr
n_pr_landing_chosen = n_pr

CL_max_clean_chosen = 1.65
CL_max_climb_chosen = 1.65
CL_max_to_chosen = 1.95
CL_max_land_chosen = 2.3
%% Performance calculations (Deprecated section. No need to check over this. Move onto detailed performance calculations)

% Harmonic Range (Take off at MPL and MTOW)
h_cruise_harmonic = 22000*0.3048
W0_harmonic = W0_MTOW_chosen
W_payload_harmonic = W_payload_max
Wf_harmonic = W0_harmonic - W_payload_harmonic - We_chosen - W_batteries_chosen
recharge_reserve_harmonic = 0.05
reserve_fuel_harmonic = Wf_harmonic*fuel_reserve
recharge_fuel_harmonic = Wf_harmonic*recharge_reserve_harmonic
Wf_available_harmonic = Wf_harmonic-reserve_fuel_harmonic-recharge_fuel_harmonic
W1_harmonic = W1_W0_constraint*W0_harmonic
Wf_TO_harmonic = W0_harmonic-W1_harmonic
W2_harmonic = W2_W1_constraint*W1_harmonic
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
W3_harmonic = (W2_harmonic-Wf_available_harmonic+Wf_TO_harmonic + Wf_climb_harmonic)/(2-W4_W3_constraint)
R_h_chosen = n_pr_cruise_chosen/(g*cp_cruise_chosen)*L_D_cruise_aim_constraint*log(W2_harmonic/W3_harmonic) %km
R_h_chosen = R_h_chosen/1000/1.852
R_h_chosen = double(subs(R_h_chosen,h,h_cruise_harmonic))

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
TOP_chosen = P0_W0_chosen^-1*W0_S_chosen*1/sigma_to*1/CL_max_to
TOFL_chosen = 11.8*TOP_chosen + 0.255*TOP_chosen^2
   
    

% Service ceiling
sigma_climb_chosen = sigma
alpha_climb_chosen = 0.9*sigma_climb_chosen^0.8
beta_climb_chosen = W1_W0*W2_W1

eq_2 = 100*0.3048/60 == n_pr*alpha_climb_chosen/beta_climb_chosen*P0_W0_chosen-2/(sigma_climb_chosen*rho_SL)*(K/(3*CD0_climb))^(1/2)*(W0_S_chosen)^(1/2)*(1.155*beta_climb_chosen^(1/2)/L_D_star)
assume(h,'real')
service_ceiling_chosen = double(vpasolve(eq_2,h))
service_ceiling_chosen = service_ceiling_chosen/0.3048

% Absolute ceiling
eq_3 = 0 == n_pr*alpha_climb_chosen/beta_climb_chosen*P0_W0_chosen-2/(sigma_climb_chosen*rho_SL)*(K/(3*CD0_climb))^(1/2)*(W0_S_chosen)^(1/2)*(1.155*beta_climb_chosen^(1/2)/L_D_star)
assume(h,'real')
absolute_ceiling = double(vpasolve(eq_3,h))
absolute_ceiling = absolute_ceiling/0.3048


%% Detailed performance calculations (Takeoff and climb)
% Setting intial weight
W0_chosen = W0_MTOW_chosen
M0_chosen = M0_MTOW_chosen
% Total fuel weight
Wf_chosen = W0_chosen - We_chosen - W_batteries_chosen - W_payload_max
Wf_reserve = fuel_reserve*Wf_chosen
Wf_available = Wf_chosen - Wf_reserve
% Height to cruise at
h_cruise_chosen = 16000*0.3048

% W and M is the instantaneous weight
W = W0_chosen
M = M0_chosen
% Instantaneous Battery energy
E_total_battery = E_total_battery_chosen
E_battery_1 = E_battery_1_chosen
E_battery_2 = E_battery_2_chosen
E_battery_3 = E_battery_3_chosen
E_battery_4 = E_battery_4_chosen


% Take off
% Power is at max
P_ICE_to = P0_ICE_chosen;
P_elec_to = P0_elec_chosen;
P_to = P_ICE_to + P_elec_to;
% For ground roll (follows Set 6 slides (Toby))
k_lof_chosen = 1.1
CL0_to = 0.1 % From Qunyh's CL-AoA graph
V_lof_chosen = k_lof_chosen*(2*W0_chosen/(rho_SL*CL_max_to_chosen*S_chosen))^(1/2)
% Friction coefficient
% u_to = 0.05 % For hard ground
% u_to = 0.1 % For short dry grass
u_to = 0.3 % For long wet grass
% Static thrust
T_static = n_pr_to_chosen*P_to/V_lof_chosen
h_wing_chosen = 2.7 %m % Height of wing above the ground
K_IGE_to = 16*(h_wing_chosen/b_chosen)^2/(1+16*(h_wing_chosen/b_chosen)^2)*K_to_chosen; % Adjusting K for ground effect
fuel_mass_flow_rate_to = cp_to_chosen*P_ICE_to;

V = 0;
s_to_g =  0;
dt_to = 0.01
iterate_limit = 3600/dt_to; % Limiting max takeoff to be an hour before loop breaks
CL_to_g = CL0_to;
time_to_g = 0;
% For the calculated static thrust, time step through the ground roll TO
for iterate = 1:iterate_limit
    V_stall_to_g_array(iterate) =  (2*W/S_chosen/(rho_SL*CL_max_to_chosen))^(1/2); % Just getting the V_stall, doesnt effect other variables
    V = V + (T_static-1/2*rho_SL*V^2*S_chosen*(CD0_to_chosen+K_IGE_to*CL_to_g^2)-u_to*(W-1/2*rho_SL*V^2*S_chosen*CL_to_g))/M*dt_to;
    s_to_g = s_to_g + V*dt_to;
    time_to_g = time_to_g + dt_to;
    W = W-fuel_mass_flow_rate_to*g*dt_to;
    E_total_battery = E_total_battery - P_elec_to*dt_to/n_e;
    E_battery_1 = E_battery_1 - P_elec_to*dt_to/n_e/2;
    E_battery_2 = E_battery_2 - P_elec_to*dt_to/n_e/2;
    M = W/g;
    if V > V_lof_chosen % Break if velocity exceeds lift off velocity
        break
    end
end

% Setting results at end of ground roll
V_stall_to_g = mean(V_stall_to_g_array)
s_to_g_calculated = s_to_g
time_to_g_calculated = time_to_g
W_end_to_g = W
M_end_to_g = M
V_end_to_g = V

E_total_battery_end_to_g = E_total_battery
E_battery_1_end_to_g = E_battery_1
E_battery_2_end_to_g = E_battery_2
E_battery_3_end_to_g = E_battery_3
E_battery_4_end_to_g = E_battery_4

Wf_used_to_g = W0_chosen - W_end_to_g
Mf_used_to_g = Wf_used_to_g/g

E_battery_used_to_g = E_total_battery_chosen - E_total_battery;


% For air distance (Set 6 slides (Toby))
CL_alpha_to = 5.53859 % From slope Quynh sent
pitch = 30 %deg
pitch = pitch*pi/180 %rad


k_2_chosen = 1.13
V2_chosen = k_2_chosen*(2*W0_chosen/(rho_SL*CL_max_to_chosen*S_chosen))^(1/2)
dV_tolerance_air = 0.1;
pitch_increment = 0.2;
h_clearance_chosen = 35*0.3048

pitch_iterate_limit = 1000
for pitch_iterate = 1:pitch_iterate_limit % Loop to see if pitch is correct
    s_to_air = 0;
    h_to_air = 0;
    time_to_air = 0;
    W = W_end_to_g;
    M = M_end_to_g;
    V = V_end_to_g;

    for iterate = 1:iterate_limit % Loop to calculate the performance parameters and see V is around V2
        V_stall_to_air_array(iterate) = (2*W/S_chosen/(rho_SL*CL_max_to_chosen))^(1/2);
        CL_to_air = W*cos(pitch)/(1/2*rho_SL*V^2*S_chosen)
        V = V+(T_static-1/2*rho_SL*V^2*S_chosen*(CD0_to_chosen+K_IGE_to*CL_to_air^2)-W*sin(pitch))/M*dt_to;
        s_to_air= s_to_air+V*cos(pitch)*dt_to;
%         AoA = (CL_to_air-CL0_to)/(CL_alpha_to);
%         Vv = V*sin(pitch-AoA);
        Vv = V*sin(pitch);
        h_to_air = h_to_air + Vv*dt_to;
        W = W-fuel_mass_flow_rate_to*g*dt_to;
        M = W/g;
        E_total_battery = E_total_battery - P_elec_to*dt_to/n_e;
        E_battery_1 = E_battery_1 - P_elec_to*dt_to/n_e/2;
        E_battery_2 = E_battery_2 - P_elec_to*dt_to/n_e/2;
        time_to_air = time_to_air + dt_to;
        if h_to_air > h_clearance_chosen
            break
        end

    end
    if V < V2_chosen && abs(V-V2_chosen) > dV_tolerance_air
        pitch = pitch - pitch_increment*pi/180;
    elseif V > V2_chosen && abs(V-V2_chosen) > dV_tolerance_air
        pitch = pitch + pitch_increment*pi/180;
    else
        break
    end
end
% Setting results at end of air TO
V_stall_to_air = mean(V_stall_to_air_array)

W_end_to_air = W
M_end_to_air = M

E_total_battery_end_to_air = E_total_battery
E_battery_1_end_to_air = E_battery_1
E_battery_2_end_to_air = E_battery_2
E_battery_3_end_to_air = E_battery_3
E_battery_4_end_to_air = E_battery_4

Wf_used_to_air = W_end_to_g - W_end_to_air
Mf_used_to_air = Wf_used_to_air/g

E_battery_used_to_air = E_total_battery_end_to_g - E_total_battery
W1_calculated = W_end_to_air
M1_calculated = M_end_to_air

s_to_air_calculated = s_to_air
s_total = s_to_g_calculated + s_to_air_calculated

% Climb (using basic climb found in Raymer performance chapter, tried to adapt Toby's notes but it lead to nonsensical results)
% Finding the best climb rate at different height blocks
CL0_climb = 0.1
CL_alpha_climb = 5.53859
% Break up climb into blocks
no_blocks = 4;
h_cruise_blocks = linspace(h_clearance_chosen,h_cruise_chosen,no_blocks+1);

alpha_climb_chosen = sigma^0.8; % From set 6 slides (Hugh)


% Tolerances
Vv_tolerance = 0.01

% Finding the best climb rate in these blocks
time_climb = 0
s_climb = 0
P_elec_climb = P0_elec_chosen;
gamma_tolerance = 0.001;
for iBlock = 1:no_blocks % Looping over different segments of climb
    % Getting climb rates for a range of velocities
    alpha_climb_block = double(subs(alpha_climb_chosen,h,mean([h_cruise_blocks(iBlock) h_cruise_blocks(iBlock + 1)])));
    sigma_climb_block = double(subs(sigma,h,mean([h_cruise_blocks(iBlock) h_cruise_blocks(iBlock + 1)])));
    rho = sigma_climb_block*rho_SL;

    V_stall = (2*W/S_chosen/(rho*CL_max_climb_chosen))^(1/2);
    V_range = linspace(V_stall,150,100);

    for iV = 1:length(V_range) % Looping over all the velocities
        iterate_limit = 1000;
        gamma_block = 0;
        for iterate = 1:iterate_limit % Loop to solve for the climb angle (commented out sections was just a slow symbolic way of solve this)
        % If thrust is in direction of the flight path
    %         syms gamma_block_sym
            V = V_range(iV);
            P_ICE_climb = alpha_climb_block*P0_ICE_chosen;
            P_climb = P_ICE_climb + P_elec_climb;
            T = alpha_climb_block*n_pr*P_climb/V;
%             L = W*cos(gamma_block_sym);
            L = W*cos(gamma_block);
            CL = L/(1/2*rho*V^2*S_chosen);
            D = 1/2*rho*V^2*S_chosen*(CD0_climb_chosen+K_climb_chosen*CL^2);
    %         eq1 = sin(gamma_block_sym) == (T-D)/W;
%             assume(gamma_block_sym > 0 & gamma_block_sym < pi/2)
%             gamma_block = vpasolve(eq1,gamma_block_sym,[0,pi/2]);
%             gamma_block = double(gamma_block);
            gamma_block_new = asin((T-D)/W);
            if abs(gamma_block_new - gamma_block) < gamma_tolerance
                break
            end
            gamma_block = gamma_block_new;
        end
        % Calculating climb rate
        Vv = V*sin(gamma_block); 
        % Storing results in an array
        AoA_block_array(iV) = (CL-CL0_climb)/CL_alpha_climb;
        pitch_block_array(iV) = AoA_block_array(iV)+gamma_block;
        gamma_block_array(iV) = gamma_block;
        Vv_block_array(iV) = Vv;
    end
    % Getting best climb rate by taking the max climb rate and the other
    % parameters with it
    [Vv_best_block, best_index] = max(Vv_block_array,[],'all') ;
    V_best(iBlock) = V_range(best_index);
    AoA_best(iBlock) = AoA_block_array(best_index);
    pitch_best(iBlock) = pitch_block_array(best_index);
    gamma_best(iBlock) = gamma_block_array(best_index);
%     L_D_best(iBlock) = L_D_block_array(best_index);
    
    % Climbing at that best speed associated with the best climb rate over the segment
    height_discretization = linspace(h_cruise_blocks(iBlock), h_cruise_blocks(iBlock+1),100); % Finer height discretization
    dh = height_discretization(2) - height_discretization(1);
    V = V_best(iBlock);
    % Calculating the climb rate at the best speed for that given height
    for iHeight = 1:length(height_discretization) 
        alpha_climb_block_recalc = double(subs(alpha_climb_chosen,h,height_discretization(iHeight)));
        sigma_climb_block_recalc = double(subs(sigma,h,height_discretization(iHeight)));
        rho = sigma_climb_block*rho_SL;
        P_ICE_climb = alpha_climb_block_recalc*P0_ICE_chosen;
        P_climb = P_ICE_climb + P_elec_climb;
        gamma_block_recalc = 0;
        for iterate = 1:iterate_limit
            T = alpha_climb_block*n_pr*P_climb/V;
    %         L = W*cos(gamma_block_sym);
            L = W*cos(gamma_block_recalc);
            CL = L/(1/2*rho*V^2*S_chosen);
            D = 1/2*rho*V^2*S_chosen*(CD0_climb_chosen+K_climb_chosen*CL^2);
%             eq1 = sin(gamma_block_sym) == (T-D)/W;

%             assume(gamma_block_sym > 0 & gamma_block_sym < pi/2)
%             gamma_block_recalc = vpasolve(eq1,gamma_block_sym,[0,pi/2]);
%             gamma_block_recalc = double(gamma_block_recalc);
            gamma_block_recalc_new = asin((T-D)/W);
            if abs(gamma_block_recalc_new - gamma_block_recalc) < gamma_tolerance
                break
            end
            gamma_block_recalc = gamma_block_recalc_new;
        end
        Vv = V*sin(gamma_block_recalc);
        % Calculate other climb paramters for this height
        AoA_block_recalc_array((iBlock-1)*length(height_discretization) + iHeight) = (CL-CL0_climb)/CL_alpha_climb;
        pitch_block_recalc_array((iBlock-1)*length(height_discretization) + iHeight) = AoA_block_recalc_array((iBlock-1)*length(height_discretization) + iHeight)+gamma_block;
        gamma_block_recalc_array((iBlock-1)*length(height_discretization) + iHeight) = gamma_block;
        V_stall_block(iHeight) = (2*W/S_chosen/(rho*CL_max_climb_chosen))^(1/2);
        % Calculate progression of climb
        dt_climb = dh/Vv;
        time_climb = dt_climb + time_climb;
        s_climb = s_climb + V*cos(gamma_block_recalc);
        fuel_mass_flow_rate_climb = alpha_climb_block*P0_ICE_chosen*cp_climb_chosen;
        W = W - fuel_mass_flow_rate_climb*g*dt_climb;
        E_total_battery = E_total_battery - P_elec_climb*dt_climb/n_e;
        E_battery_1 = E_battery_1 - P_elec_climb*dt_climb/n_e/2;
        E_battery_2 = E_battery_2 - P_elec_climb*dt_climb/n_e/2;
    end
    V_stall_climb(iBlock) = mean(V_stall_block);
end
% Setting results after climb
W_end_climb = W;

W2_calculated = W_end_climb
M2_calculated = W_end_climb/g

Wf_used_climb = W1_calculated - W2_calculated
Mf_used_climb = Wf_used_climb/g

E_total_battery_end_climb = E_total_battery;
E_battery_1_end_climb = E_battery_1
E_battery_2_end_climb = E_battery_2
E_battery_3_end_climb = E_battery_3
E_battery_4_end_climb = E_battery_4

E_battery_used_climb =  E_total_battery_end_to_air - E_total_battery


%% Detailed performance calculations (Cruise)
% Cruise conditions
V_cruise_chosen = V_cruise_req
dt_cruise = 0.1
sigma_cruise_chosen = double(subs(sigma,h,h_cruise_chosen))
alpha_cruise_chosen = 0.7*sigma_cruise_chosen^0.8
rho_cruise_chosen = sigma_cruise_chosen*rho_SL

% Percentage levels of the two batteries at the start of cruise
E_total_battery_start_cruise = E_total_battery_end_climb;
E_battery_1_start_cruise = E_battery_1_end_climb;
E_battery_2_start_cruise = E_battery_2_end_climb;
E_battery_1_start_cruise_percentage = E_battery_1_start_cruise/E_battery_1_chosen;
E_battery_2_start_cruise_percentage = E_battery_2_start_cruise/E_battery_2_chosen;

E_total_battery = E_total_battery_start_cruise;
E_battery_1 = E_battery_1_start_cruise;
E_battery_2 = E_battery_2_start_cruise;
E_battery_3 = E_battery_3_chosen;
E_battery_4 = E_battery_4_chosen;

% Power available
P_available_ICE = P0_ICE_chosen*alpha_cruise_chosen
% P_available_elec = 0.02*P0_ICE_chosen
% Power available coming out of the motors
P_available_elec = 22000
% Power recharing going into the battery
P_recharge = 20000 % W
% Total power available
P_available = P_available_ICE + P_available_elec - P_recharge/n_e;

iterate_limit = 10000000
s_cruise = 0;
time_cruise = 0;
W = W2_calculated;

% Plotting power required to power available chart just to see if the
% height I'm working at is reasonable (can disable)
plot_power_chart = true;
if plot_power_chart == true
    alpha_cruise_chart = 0.7*sigma^0.8;
    P_available_ICE_chart = P0_ICE_chosen*alpha_cruise_chart;
    P_available_chart = P_available_ICE_chart + P_available_elec - P_recharge/n_e;
    rho_chart = rho_SL*sigma;
    CL_chart = W/(1/2*rho_chart*V_cruise_chosen^2*S_chosen);
    T_chart = 1/2*rho_chart*V_cruise_chosen^2*S_chosen*(CD0_cruise_chosen+K_cruise_chosen*CL_chart^2);
    P_cruise_chart = T_chart*V_cruise_chosen/n_pr_cruise_chosen;
    h_axis = 10:100:25000*0.3048;
    figure
    plot(h_axis/0.3048, double(subs(P_available_chart,h,h_axis)))
    hold on
    plot(h_axis/0.3048, double(subs(P_cruise_chart,h,h_axis)))
    xlabel('h (ft)')
    ylabel('Power (W)')

    legend('Power available','Power required')
end

% Recharging states (set 1 is battery 1 and 2 and set 2 is battery 3 and 4)
recharging_set_1 = true;
recharging_set_2 = false;
% Percent to stop recharging
percent_to_recharge = 1;

switch_iterate_set_1 = 1;
switch_iterate_set_2 = 1;
% Time when a switch happens
time_finished_recharging_set_1 = [];
time_finished_recharging_set_2 = [];
% Timestep through cruise
for iterate = 1:iterate_limit
    V_stall_cruise_array(iterate) = (2*W/S_chosen/(rho_cruise_chosen*CL_max_clean_chosen))^(1/2);
    CL = W/(1/2*rho_cruise_chosen*V_cruise_chosen^2*S_chosen);
    T = 1/2*rho_cruise_chosen*V_cruise_chosen^2*S_chosen*(CD0_cruise_chosen+K_cruise_chosen*CL^2);
    % Power required
    P_cruise = T*V_cruise_chosen/n_pr_cruise_chosen;
    % Electrical power at cruise
    P_elec_cruise = P_available_elec;
    if P_elec_cruise < P_recharge
        disp('Not worth using electric at this recharge rate')
        break
    end
    % Power from combustion engine to fly
    P_ICE_to_fly = P_cruise - P_available_elec;
    if P_available < P_cruise
        disp('Unable to cruise at this altitude')
        break
    end
    % Excess power for optimization (decreases throughout cruise as W
    % decreases)
    P_excess(iterate) = P_available - P_cruise;
    fuel_mass_flow_rate_cruise = cp_cruise_chosen*(P_ICE_to_fly + P_recharge/n_e);
    W = W-fuel_mass_flow_rate_cruise*dt_cruise*g;

    s_cruise = s_cruise + dt_cruise*V_cruise_chosen;
    time_cruise = time_cruise + dt_cruise;
    
    % Recharge cycle
    if recharging_set_1 == true
        E_battery_1 = E_battery_1 + P_recharge/no_batteries/2*dt_cruise;
        E_battery_2 = E_battery_2 + P_recharge/no_batteries/2*dt_cruise;
        E_battery_3 = E_battery_3 - P_elec_cruise/n_e*dt_cruise;
        E_battery_4 = E_battery_4 - P_elec_cruise/n_e*dt_cruise;
    elseif recharging_set_2 == true
        E_battery_3 = E_battery_3 + P_recharge/no_batteries/2*dt_cruise;
        E_battery_4 = E_battery_4 + P_recharge/no_batteries/2*dt_cruise;
        E_battery_1 = E_battery_1 - P_elec_cruise/n_e*dt_cruise;
        E_battery_2 = E_battery_2 - P_elec_cruise/n_e*dt_cruise;
    end
    
    % Recharge cycle switching condition
    if E_battery_1 >= percent_to_recharge*E_battery_1_chosen && E_battery_2 >= percent_to_recharge*E_battery_2_chosen
        recharging_set_2 = true;
        recharging_set_1 = false;
        time_finished_recharging_set_1(switch_iterate_set_1) = time_cruise;
        switch_iterate_set_1 = switch_iterate_set_1 + 1;
    elseif E_battery_3 >= percent_to_recharge*E_battery_3_chosen && E_battery_4 >= percent_to_recharge*E_battery_4_chosen
        recharging_set_2 = false;
        recharging_set_1 = true;
        time_finished_recharging_set_2(switch_iterate_set_2) = time_cruise;
        switch_iterate_set_2 = switch_iterate_set_2 + 1;
    end

    if E_battery_1 < 0 || E_battery_2 < 0 || E_battery_3 < 0 || E_battery_4 < 0
        disp('Parameters lead to complete depletion of batteries')
        break
    end

    if W < W0_chosen-Wf_available % Stop when we have used all our fuel available
        break
    end
end

% Setting results at end of cruise
V_stall_cruise = mean(V_stall_cruise_array);

s_cruise_calculated = s_cruise;
s_cruise_calculated = s_cruise/1000/1.852; % nm
W_end_cruise = W;

W3_calculated = W;
M3_calculated = W3_calculated/g;

Wf_used_cruise = W2_calculated - W3_calculated;
Mf_used_cruise = Wf_used_cruise/g;

disp('At end of flight:')
disp(num2str(E_battery_1/E_battery_1_chosen*100,'Battery 1 has %.2f percent left'))
disp(num2str(E_battery_2/E_battery_2_chosen*100,'Battery 2 has %.2f percent left'))
disp(num2str(E_battery_3/E_battery_3_chosen*100,'Battery 3 has %.2f percent left'))
disp(num2str(E_battery_4/E_battery_4_chosen*100,'Battery 4 has %.2f percent left'))
disp(num2str(P_excess(1)/1000,'Initial excess power at this condition is %.2f kW'))

%% Detailed performance calculations (Descent, Approach and Landing)
% Descent
V_descent = V_cruise_chosen
gamma_descent = 3
gamma_descent = gamma_descent*pi/180
h_land = 50 %ft
h_land = h_land*0.3048 %m
W = W3_calculated
Vv_descent = V_descent*sin(gamma_descent)
h_descent = h_cruise_chosen
s_descent = 0
t_descent = 0
dt_descent = 1
iterate_limit = 3600/dt_descent
for iterate = 1:iterate_limit
    rho = double(subs(sigma,h,h_descent));
    L_descent = W*cos(gamma_descent);
    CL_descent = 2*L/(rho*V_descent^2);
    CD = CD0_descent_chosen + K_descent_chosen*CL^2;
    D = 1/2*rho*V_descent^2*CD;
    T = -W*sin(gamma_descent)+D;
    P = T*V_descent/n_pr_descent_chosen;
    fuel_mass_flow_rate_descent = cp_descent_chosen*P;
    W = W - fuel_mass_flow_rate_descent*g*dt_descent;
    t_descent = t_descent + dt_descent;
    h_descent = h_descent- dt_descent*Vv_descent;
    s_descent = s_descent + dt_descent*V_descent*cos(gamma_descent);
    if h_descent < h_land
        break
    end
end
W_end_of_descent = W
W4_calculated = W_end_of_descent
Wf_used_descent = W4_calculated - W3_calculated % Increases with increasing gamma_descent

% Landing

% Approach/Air
W = W4_calculated
gamma_approach = 3
gamma_approach = gamma_approach*pi/180
k_app_chosen = k_land
V_stall_land = (2/rho_SL*W/S_chosen*1/CL_max_land_chosen)^(1/2)
V_approach = k_app_chosen*V_stall_land
s_land_air = h_land/tan(gamma_approach)

% Ground roll
s_land_g = 0
CL0_land = 0.1
CL = CL0_land
k_td = 1.15
V_td = k_td*V_stall_approach
V = V_td
dt_land_g = 0.01
iterate_limit = 3600/dt_descent
K_IGE_land = 16*(h_wing_chosen/b_chosen)^2/(1+16*(h_wing_chosen/b_chosen)^2)*K_landing_chosen; % Adjusting K for ground effect
CD = CD0_landing_chosen+K_IGE_land*CL^2
u_land = 0.3
s_land_g = 0
t_land_g = 0
for iterate = 1:iterate_limit
    L = 1/2*rho_SL*V^2*S_chosen*CL;
    D = 1/2*rho_SL*V^2*S_chosen*CD;
    N = W-L;
    Fr = u_land*N;
    M = W/g;
    V = V +(-D-Fr)/M*dt_land_g;
    s_land_g = s_land_g + V*dt_land_g;
    t_land_g = t_land_g + dt_land_g;
    if V < 0
        break
    end

end
s_land = s_land_g + s_land_air




%% Other parameter calculations
% Service ceiling
sigma_SC_chosen = sigma
alpha_SC_chosen = 0.9*sigma^0.8
W_SC = W0_MTOW_chosen

% If all engines running
% P_SC_ICE = alpha_SC_chosen*P0_ICE_chosen
% P_SC_elec = P0_elec_chosen

% If all motors fail (maybe due to an electrical short)
% P_SC_ICE = alpha_SC_chosen*P0_ICE_chosen
% P_SC_elec = 0

% If both ICE flame out
% P_SC_ICE = 0
% P_SC_elec = P0_elec_chosen

% If one engine fails completely
P_SC_ICE = alpha_SC_chosen*P0_ICE_chosen/2
P_SC_elec = P0_elec_chosen/2


P_SC = P_SC_ICE + P_SC_elec
rho_SC = sigma_SC_chosen*rho_SL
L_D_star_SC = 1/(4*CD0_climb_chosen*K_climb_chosen)^(1/2)

eq_2 = 100*0.3048/60 == n_pr_climb_chosen*P_SC/W_SC-2/(rho_SC)*(K_climb_chosen/(3*CD0_climb_chosen))^(1/2)*(W_SC/S_chosen)^(1/2)*(1.155/L_D_star_SC)
assume(h,'real')
service_ceiling_chosen = double(vpasolve(eq_2,h))
service_ceiling_chosen = service_ceiling_chosen/0.3048



%% Summarizing results (old code)
summarize_results = false
if summarize_results == true
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
    disp(num2str(Mf_chosen_harmonic,'Fuel mass = %.2f kg (for harmonic range)'))
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
end







