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