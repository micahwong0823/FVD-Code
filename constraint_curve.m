% Non-constraint analysis (because Jayson says that constraint analysis is overkill but can be used for the appendices)
% For takeoff
% Figure formatting
set(0,'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize',10,...
    'DefaultAxesFontSize',16,...
    'DefaultTextInterpreter','latex',...
    'DefaultLegendInterpreter','latex',...
    'DefaultLegendFontSize',16)

% Using estimates of W0_S from Raymer
W0_S_arr = linspace(2500,3500,30)

m_batt_arr = []
P0_W0_cr_arr = []
P0_W0_sc_arr = []
P0_W0_to_arr = []
R_cr_nm = []
% P0_elec_P0_opt_arr = []
h_cr_arr = []

for iW0_S = 1:length(W0_S_arr)
    W0_S = W0_S_arr(iW0_S)

    % Using preliminary estimated mass
    M0 = 1.949591926490346e+04 %kg
    %Gravity
    g = 9.8065
    % Weight
    W0 = M0*g

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

    % Draft calculation on time taken (wildly off)
%     % Takeoff factor on stall speed
%     k_to = 1.2
%     % Lift at takeoff
%     Cl_2 = Cl_max_to/k_to^2
%     % Sea level density
%     rho_SL = 1.225
%     % Velocity at takeoff
%     V_2 = (2/rho_SL*W0_S*1/Cl_2)^(1/2) %m/s
%     % Propeller blades
%     K_p = 0.49
%     % Blade diameter and area
%     D_p = K_p*(P0_to/1000)^(1/4)
%     A_p = (pi*D_p^2/4)
%     % Vary velocity
%     syms V
%     % First segment of thrust
%     T_seg_1 = (2*P0_to^2*rho_SL*A_p)^(1/3)
%     % Propulsive efficiency
%     n_p = 0.85
%     % Second segment of thrust
%     T_seg_2 = (n_p*n_e*P0_elec_to+n_p*P0_ICE_to)/V
%     % Intersecting velocity between T_seg_1 and 2
%     V_int = double(solve(T_seg_1==T_seg_2,V))
%     % Integrate the thrust over V = 0 to V = V2
%     T_integrated = double(int(T_seg_1,V,[0,V_int]) + int(T_seg_2,V,[V_int,V_2]))
%     % Average acceleration
%     a_average = T_integrated/M0*1/V_2
%     % Reduce
%     r_t = 0.8
%     a_average = 0.8*a_average
%     % Roll runway
%     s_run = V_2^2/(2*a_average)
%     % Clearance height
%     h_to = 35 %ft
%     h_to = h_to*0.3048
%     T_2 = double(subs(T_seg_2,V,V_2))
%     % Aerodynamic parameters
%     Cd0_to = 0.02 % Torenbeek
%     e = 0.7
%     % Aspect Ratio
%     A = 12
%     K = 1/(pi*A*e)
%     % Drag at end of takeoff
%     Cd_2 = Cd0_to+K*Cl_2^2
%     % Weight factor
%     beta_to = 0.97
%     % Takeoff climb angle
%     tan_gamma_2 = (T_2/(beta_to*W0)-Cd_2/Cl_2)
%     gamma_2 = atan(tan_gamma_2)*180/pi
%     s_air = h_to/tan_gamma_2
%     % Takeoff
%     s_to = s_run + s_air
%     % Adding 15% margin
%     s_to = 1.15*s_to
% 
%     % Time taken
%     t_run = V_2/a_average
%     t_air = h_to/(V_2*sin(gamma_2*pi/180))
%     t_to = t_run + t_air
%     % Voltage of motor
%     V_motor = 750 %V
%     % Battery capacity used
%     E_to = t_to*P0_elec_to
%     % Battery capacity left
%     E_left = E_batt-E_to


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
    % Efficiency factor
    e = 0.7
    % Aspect Ratio
    A = 12
    K = 1/(pi*A*e)
    % Max lift to drag ratio
    L_D_star = 1/(4*Cd0_cr*K)^(1/2)
    % Density sea level
    rho_SL = 1.225
    % Lift to drag ratio
    Cl_cr = 2/(sigma*rho_SL*V_cr_req^2)*W0_S
    Cd_cr = Cd0_cr + K*Cl_cr^2
    L_D_cr = Cl_cr/Cd_cr
    % PSFC
    cp_cr = 0.085*10^-6 % kg/Ws
    cp_cr = 0.7*cp_cr % Reduce due to parallel
    % Gravity
    g = 9.8065

    % Propulsive efficiency
    n_p = 0.85
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
%     figure
%     fplot(P0_cr)
%     xlim([0 25000*0.3048])


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

%     % If cruise is run on fuel then calculate the percentage of electricity
%     % used for takeoff
%     P0_ICE_P0_opt = P0_cr_opt/P0_to
%     P0_elec_P0_opt = 1-P0_ICE_P0_opt

    % Power distribution
    % Power motor for to
    P0_elec_to = P0_elec_P0*P0_to
    % Power ICE for to
    P0_ICE_to = P0_ICE_P0*P0_to

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


    % Range
    sigma_cr_opt = double(subs(sigma,h,h_cr))
    W3_W2_cr_opt = double(subs(W3_W2,h,h_cr))
    beta_cr_opt = subs(beta_cr,h,h_cr)
    Cl_cr_opt = 2/(rho_SL*sigma_cr_opt*V_cr_req^2)*beta_cr_opt*W0_S
    Cd_cr_opt = Cd0_cr + K*Cl_cr_opt^2
    L_D_cr_opt = Cl_cr_opt/Cd_cr_opt
    R_cr = double(n_p/(g*cp_cr)*L_D_cr_opt*log(W3_W2_cr_opt^-1)) %m
    R_cr_nm = R_cr/1.852/1000

    % Getting constraint curve and other important information
    P0_W0_to = P0_to/W0
    P0_W0_cr = P0_cr_opt/W0
    P0_W0_sc = P0_sc/W0

    
    m_batt_arr(iW0_S) = m_batt
    P0_W0_cr_arr(iW0_S) = P0_W0_cr
    P0_W0_sc_arr(iW0_S) = P0_W0_sc
    P0_W0_to_arr(iW0_S) = P0_W0_to
    R_cr_nm(iW0_S) = R_cr_nm
%     P0_elec_P0_opt_arr(iW0_S) = P0_elec_P0_opt
    h_cr_arr(iW0_S) = h_cr

end
constraint_plt = figure
plot(W0_S_arr,P0_W0_to_arr)
hold on
plot(W0_S_arr,P0_W0_cr_arr)
hold on
plot(W0_S_arr,P0_W0_sc_arr)
xlabel('$\frac{W_0}{S}$')
ylabel('$\frac{P_0}{W_0}$')
title('Constraint Curve (FP)')
legend('Takeoff','Cruise','Service Ceiling')
grid on


m_batt_plt = figure
plot(W0_S_arr,m_batt_arr)
xlabel('$\frac{W_0}{S}$')
ylabel('$m_{batt}$')
title('Battery mass curve')
grid on
% ylim([0 Inf])

% elec_perc_plt = figure
% plot(W0_S_arr,P0_elec_P0_opt_arr)
% xlabel('$\frac{W_0}{S}$')
% ylabel('$\frac{P_{0_e}}{P_0}$')
% title('Electricity fraction')
% grid on
% ylim([0 Inf])
% % R_cr_nm_plt = figure
% % plot(W0_S_arr,R_cr_nm)

h_cr_plt = figure
plot(W0_S_arr,h_cr_arr)
xlabel('$\frac{W_0}{S}$')
ylabel('h (m)')
title('Optimal cruise height')
grid on

