% Old design code

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

%% Old detailed climb attempt
%         AoA_block = 0.01;
%         pitch_block = 0.5;
%         gamma_block = 0.01;
%         Vv = 0;
%         for iterate = 1:iterate_limit
%             V = V_range(iV);
%             T = alpha_climb_block*n_pr*P0_chosen/V;
% %             T = 0
%             L = W*cos(gamma_block)-T*sin(AoA_block);
%             CL = L/(1/2*rho*V^2*S_chosen);
%             D = 1/2*rho*V^2*S_chosen*(CD0+K*CL^2);
% %             gamma_block = atan((T*cos(pitch_block)*sin(gamma_block)-L+W*cos(gamma_block))/(T*cos(AoA_block)-D - W*sin(gamma_block)));
%             gamma_block = atan((T*sin(pitch_block)-W+L*cos(gamma_block)-D*sin(gamma_block))/(T*cos(pitch_block)-L*sin(gamma_block)-D*cos(gamma_block)));
%             AoA_block = (CL-CL0_climb)/CL_alpha_climb;
%             pitch_block = gamma_block + AoA_block;
%             Vv_new = V*sin(gamma_block);
%             if abs(Vv_new - Vv) < Vv_tolerance
%                 break
%             end
%             Vv = Vv_new;
%         end

%         syms AoA_block_sym pitch_block_sym gamma_block_sym
%         V = V_range(iV);
%         T = alpha_climb_block*n_pr*P0_chosen/V;
%         L = W*cos(gamma_block_sym)-T*sin(AoA_block_sym);
%         CL = L/(1/2*rho*V^2*S_chosen);
%         D = 1/2*rho*V^2*S_chosen*(CD0+K_climb_chosen*CL^2);
%         eq1 = AoA_block_sym == (CL-CL0_climb)/CL_alpha_climb;
%         eq2 = pitch_block_sym == gamma_block_sym + AoA_block_sym;
%         eq3 = tan(gamma_block_sym) == (T*sin(pitch_block_sym)-W+L*cos(gamma_block_sym)-D*sin(gamma_block_sym))/(T*cos(pitch_block_sym)-L*sin(gamma_block_sym)-D*cos(gamma_block_sym));
%         
%         assume(AoA_block_sym > 0 & AoA_block_sym < pi/2)
%         assume(pitch_block_sym > 0 & pitch_block_sym < pi/2)
%         assume(gamma_block_sym > 0 & gamma_block_sym < pi/2)
% 
%         [AoA_block pitch_block gamma_block] = vpasolve([eq1 eq2 eq3],[AoA_block_sym, pitch_block_sym, gamma_block_sym]);
%         AoA_block = double(AoA_block);
%         pitch_block = double(pitch_block);
%         gamma_block = double(gamma_block);

%         syms AoA_block_sym pitch_block_sym gamma_block_sym
%         V = V_range(iV);
% %         T = alpha_climb_block*n_pr*P0_chosen/V;
% %         L = W*cos(gamma_block_sym)-T*sin(AoA_block_sym);
%         CL = W*cos(pitch_block_sym)/(1/2*rho*V^2*S_chosen);
% %         D = 1/2*rho*V^2*S_chosen*(CD0+K_climb_chosen*CL^2);
%         Vv = (n_pr_climb_chosen*alpha_climb_block*P0_chosen-(CD0-K_climb_chosen*CL^2)*V)/W;
%         
%         eq1 = AoA_block_sym == (CL-CL0_climb)/CL_alpha_climb;
%         eq2 = pitch_block_sym == gamma_block_sym + AoA_block_sym;
%         eq3 = sin(gamma_block_sym) == Vv/V;
%         
%         assume(AoA_block_sym > 0 & AoA_block_sym < pi/2)
%         assume(pitch_block_sym > 0 & pitch_block_sym < pi/2)
%         assume(gamma_block_sym > 0 & gamma_block_sym < pi/2)
% 
%         [AoA_block pitch_block gamma_block] = vpasolve([eq1 eq2 eq3],[AoA_block_sym, pitch_block_sym, gamma_block_sym]);
%         AoA_block = double(AoA_block);
%         pitch_block = double(pitch_block);
%         gamma_block = double(gamma_block);
% 
% 
%         Vv = V*sin(gamma_block);
%         AoA_block_array(iV) = AoA_block;
%         pitch_block_array(iV) = pitch_block;
%         gamma_block_array(iV) = gamma_block;
% %         L_D_block_array(iV) = L/D;
%         Vv_block_array(iV) = Vv;

%         syms AoA_block_sym pitch_block_sym gamma_block_sym
%         V = V_range(iV);
%         T = alpha_climb_block*n_pr*P0_chosen/V;
%         L = W*cos(gamma_block_sym)-T*sin(AoA_block_sym);
%         CL = L/(1/2*rho*V^2*S_chosen);
%         D = 1/2*rho*V^2*S_chosen*(CD0+K_climb_chosen*CL^2);
%         eq1 = AoA_block_sym == (CL-CL0_climb)/CL_alpha_climb;
%         eq2 = pitch_block_sym == gamma_block_sym + AoA_block_sym;
%         eq3 = tan(gamma_block_sym) == (T*cos(AoA_block_sym)-D)/(L+T*sin(AoA_block_sym));
%         
%         assume(AoA_block_sym > 0 & AoA_block_sym < pi/2)
%         assume(pitch_block_sym > 0 & pitch_block_sym < pi/2)
%         assume(gamma_block_sym > 0 & gamma_block_sym < pi/2)
% 
%         [AoA_block pitch_block gamma_block] = vpasolve([eq1 eq2 eq3],[AoA_block_sym, pitch_block_sym, gamma_block_sym]);
%         AoA_block = double(AoA_block);
%         pitch_block = double(pitch_block);
%         gamma_block = double(gamma_block);

%         AoA_block_recalc = 0;
%         pitch_block_recalc = 0;
%         gamma_block_recalc = 0; 
%         Vv = 0;
%         for iterate = 1:iterate_limit
%             T = alpha_climb_block_recalc*n_pr_climb_chosen*P0_chosen/V;
%             %             T = 0
%             L = W*cos(gamma_block_recalc)-T*sin(AoA_block_recalc);
%             CL = L/(1/2*rho*V^2*S_chosen);
%             D = 1/2*rho*V^2*S_chosen*(CD0_climb_chosen+K_climb_chosen*CL^2);
%             gamma_block_recalc = atan((T*cos(pitch_block_recalc)*sin(gamma_block_recalc)-L+W*cos(gamma_block_recalc))/(T*cos(AoA_block_recalc)-D - W*sin(gamma_block_recalc)));
%             AoA_block_recalc = (CL-CL0_climb)/CL_alpha_climb;
%             pitch_block_recalc = gamma_block_recalc + AoA_block_recalc;
%             Vv_new = V*sin(gamma_block_recalc);
%             if abs(Vv_new - Vv) < Vv_tolerance
%                 break
%             end
%             Vv = Vv_new;
%             
%         end

%% Old recharge rate code
%     if recharge_rate*time_cruise > (E_battery_used_climb + E_battery_used_to_g + E_battery_used_to_air) && abs((E_battery_used_climb + E_battery_used_to_g + E_battery_used_to_air)-recharge_rate*time_cruise) > recharge_rate_tolerance
%         recharge_rate = recharge_rate - recharge_rate_increment;
%     elseif recharge_rate*time_cruise < (E_battery_used_climb + E_battery_used_to_g + E_battery_used_to_air) && abs((E_battery_used_climb + E_battery_used_to_g + E_battery_used_to_air)-recharge_rate*time_cruise) > recharge_rate_tolerance
%         recharge_rate = recharge_rate + recharge_rate_increment;
%     else
%         break
%     end