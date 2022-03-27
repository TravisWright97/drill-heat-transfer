%ME 3371-001 Heat Transfer Project
%Travis W.

%Background info (user inputs):
L_drill = 10 * 10^-2;  %[m]
d_drill = 5 * 10^-3;  %[m]
cut_angle = 20; %[degrees]
r = d_drill / 2; %[m]
L_cone = tand(cut_angle) * (r); %[m]
h_air = 30; %[W/m^2*K], convective heat transfer coefficient of air
k_SS = 15; %[W/m*K], conductive heat transfer coefficient of stainless steel
rho_SS = 7750; %[kg/m^3], density of stainless steel
C_p_SS = 468; %[J/kg*K], specific heat capacity of stainless steel
T_inf = 20; %[C], ambient air temp/initial drill temp in Celsius

%Differential height of nodes:
del_z_1 = (L_drill - L_cone) / 6; %[m], differential length of node in z
del_z_2 = L_cone / 3; %[m], differential length of node in cone

%Generalized area and volume formulas for each kind of node
Area_1 = (1/2) * r^2; %[m^2], bottom of top node
Area_2 = r * (del_z_1/2); %[m^2], outside of top node
Area_3 = r * del_z_1; %[m^2], outside of middle-cylinder nodes
Area_4 = (1/2) * ((3*del_z_2)/(2*tand(cut_angle)))^2; %[m^2], bottom of cylinder-cone transition node
Area_5 = (1/2) * (((2*del_z_2)/(sind(cut_angle)))-((3*del_z_2)/(2*sind(cut_angle)))) * (r+((3*del_z_2)/(2*tand(cut_angle)))); %[m^2], bone contact area for cylinder-cone transition node
Area_6 = (1/2) * (del_z_2/(2*tand(cut_angle)))^2; %[m^2], bottom of middle-cone nodes
Area_7 = (1/2) * (((3*del_z_2)/(2*sind(cut_angle)))-((del_z_2)/(2*sind(cut_angle)))) * (((3*del_z_2)/(2*tand(cut_angle)))+((del_z_2)/(2*tand(cut_angle)))); %[m^2], bone contact area for middle-cone nodes
Area_8 = (1/2) * ((del_z_2)/(2*tand(cut_angle))) * ((del_z_2)/(2*sind(cut_angle))); %[m^2], bone contact area for cone-tip node
Vol_1 = (1/2) * r^2 * (del_z_1/2); %[m^3], volume of top node
Vol_2 = (1/2) * r^2 * del_z_1; %[m^3], volume of moddle-cylinder nodes
Vol_3 = ((1/2) * r^2 * (del_z_1/2)) + (((1/6)*r^2*(del_z_2/2))-((1/6)*((3*del_z_2)/(2*tand(cut_angle)))^2*((2*del_z_2)-(del_z_2/2)))); %[m^3], volume of cylinder-cone transition node
Vol_4 = ((1/6)*((3*del_z_2)/(2*tand(cut_angle)))^2*del_z_2) - ((1/6)*(del_z_2/(2*tand(cut_angle)))^2*(((3*del_z_2)/2)-del_z_2)); %[m^3], volume of middle-cone nodes
Vol_5 = (1/6) * ((del_z_2)/(2*tand(cut_angle)))^2 * (del_z_2/2); %[m^3], volume of cone-tip node

%Heat generation at cutting tip:
F_drill = 10; %[N]
Area_cone = (pi * r^2) + (pi * r * sqrt(L_cone^2 + r^2)); %[m^2]
omega_RPM = 3000; %[RPM]
omega = omega_RPM * ((2*pi)/60); %[rad/s]
Vel = r * omega; %[m/s]
mu = 0.75; 
P_drill = F_drill/Area_cone; %[Pa] 
q_fric = mu*P_drill*Vel; %[W/m^2], calculating total frictional heat gen
q_shear_1 = 0 * q_fric; %[W/m^2], diff. assumptions to find shear heat gen
q_shear_2 = 0.5 * q_fric; %[W/m^2]
q_shear_3 = 1 * q_fric; %[W/m^2]
q_gen_1 = q_fric + q_shear_1; %[W/m^2], total heat gen with different shear value assumptions
q_gen_2 = q_fric + q_shear_2; %[W/m^2]
q_gen_3 = q_fric + q_shear_3; %[W/m^2]

%Initial heat equation matrix:
T_initial = zeros(11,1) + T_inf; %initial temperatures of room temp, with extra term for T_infinity

%Time steps and errors:
del_t = 5; %[s], time step
p_1 = 0; %initial step number for 3 cases
p_2 = 0;
p_3 = 0;
t_1 = 0; %initial total time for 3 cases
t_2 = 0;
t_3 = 0;

ea_1 = 10; %initial error to get while loop started
ea_2 = 10;
ea_3 = 10;

es = 10^-6; %stopping error

T_old_1 = T_initial; %setting first temp distribution matrix to room temp for all 3 cases
T_old_2 = T_initial;
T_old_3 = T_initial;

%f and g simplification assignments:
f_1 = (h_air*Area_2*del_t)/(rho_SS*Vol_1*C_p_SS);
f_2 = (h_air*Area_2*del_t)/(rho_SS*Vol_2*C_p_SS);
f_3 = (h_air*Area_2*del_t)/(rho_SS*Vol_3*C_p_SS);

g_1 = (k_SS*Area_1*del_t)/(rho_SS*Vol_1*C_p_SS*del_z_1);
g_2 = (k_SS*Area_1*del_t)/(rho_SS*Vol_2*C_p_SS*del_z_1);
g_3_1 = (k_SS*Area_4*del_t)/(rho_SS*Vol_3*C_p_SS*del_z_2);
g_3_2 = (k_SS*Area_1*del_t)/(rho_SS*Vol_3*C_p_SS*del_z_1);
g_4_1 = (k_SS*Area_6*del_t)/(rho_SS*Vol_4*C_p_SS*del_z_2);
g_4_2 = (k_SS*Area_4*del_t)/(rho_SS*Vol_4*C_p_SS*del_z_2);
g_5 = (k_SS*Area_6*del_t)/(rho_SS*Vol_5*C_p_SS*del_z_2);

h_1_1 = (del_t*q_gen_1*Area_5)/(rho_SS*Vol_3*C_p_SS); %node 3
h_2_1 = (del_t*q_gen_1*Area_7)/(rho_SS*Vol_4*C_p_SS); %node 2
h_3_1 = (del_t*q_gen_1*Area_8)/(rho_SS*Vol_5*C_p_SS); %node 1

h_1_2 = (del_t*q_gen_2*Area_5)/(rho_SS*Vol_3*C_p_SS); %node 3
h_2_2 = (del_t*q_gen_2*Area_7)/(rho_SS*Vol_4*C_p_SS); %node 2
h_3_2 = (del_t*q_gen_2*Area_8)/(rho_SS*Vol_5*C_p_SS); %node 1

h_1_3 = (del_t*q_gen_3*Area_5)/(rho_SS*Vol_3*C_p_SS); %node 3
h_2_3 = (del_t*q_gen_3*Area_7)/(rho_SS*Vol_4*C_p_SS); %node 2
h_3_3 = (del_t*q_gen_3*Area_8)/(rho_SS*Vol_5*C_p_SS); %node 1

%A-matrix assignments:
A_1_1 = [1,0,0,0,0,0,0,0]; %assigns all the rows for A, rows/columns seperated for easier viewing/debugging of A matrices
A_1_2 = [0,g_5+1,-g_5,0,0,0,0,0];
A_1_3 = [0,-g_4_1,g_4_2+g_4_1+1,-g_4_2,0,0,0,0];
A_1_4 = [-f_3,0,-g_3_1,g_3_1+g_3_2+f_3+1,-g_3_2,0,0,0];
A_1_5 = [-f_2,0,0,-g_2,2*g_2+f_2+1,-g_2,0,0];
A_1_6 = [-f_2,0,0,0,-g_2,2*g_2+f_2+1,-g_2,0];
A_1_7 = [-f_2,0,0,0,0,-g_2,2*g_2+f_2+1,-g_2];
A_1_8 = [-f_1,0,0,0,0,0,-g_1,g_1+f_1+1];
A_1 = [A_1_1;A_1_2;A_1_3;A_1_4;A_1_5;A_1_6;A_1_7;A_1_8]; %compiles all rows into 8x8 A matrix 
A_2_1 = [1,0,0,0,0,0,0,0];
A_2_2 = [0,g_5+1,-g_5,0,0,0,0,0];
A_2_3 = [0,-g_4_1,g_4_2+g_4_1+1,-g_4_2,0,0,0,0];
A_2_4 = [-f_3,0,-g_3_1,g_3_1+g_3_2+f_3+1,-g_3_2,0,0,0];
A_2_5 = [-f_2,0,0,-g_2,2*g_2+f_2+1,-g_2,0,0];
A_2_6 = [-f_2,0,0,0,-g_2,2*g_2+f_2+1,-g_2,0];
A_2_7 = [-f_2,0,0,0,0,-g_2,2*g_2+f_2+1,-g_2];
A_2_8 = [-f_1,0,0,0,0,0,-g_1,g_1+f_1+1];
A_2 = [A_2_1;A_2_2;A_2_3;A_2_4;A_2_5;A_2_6;A_2_7;A_2_8];
A_3_1 = [1,0,0,0,0,0,0,0];
A_3_2 = [0,g_5+1,-g_5,0,0,0,0,0];
A_3_3 = [0,-g_4_1,g_4_2+g_4_1+1,-g_4_2,0,0,0,0];
A_3_4 = [-f_3,0,-g_3_1,g_3_1+g_3_2+f_3+1,-g_3_2,0,0,0];
A_3_5 = [-f_2,0,0,-g_2,2*g_2+f_2+1,-g_2,0,0];
A_3_6 = [-f_2,0,0,0,-g_2,2*g_2+f_2+1,-g_2,0];
A_3_7 = [-f_2,0,0,0,0,-g_2,2*g_2+f_2+1,-g_2];
A_3_8 = [-f_1,0,0,0,0,0,-g_1,g_1+f_1+1];
A_3 = [A_3_1;A_3_2;A_3_3;A_3_4;A_3_5;A_3_6;A_3_7;A_3_8];

i = 1; %initializing counters for the temp tracking plot 
j = 1;
k = 1;

%==========================================================================
%CASE 1 (q_shear = 0 * q_fric):

while ea_1 > es
    p_1 = p_1 + 1; %time step counter
    t_1 = p_1 * del_t; %current time value
    c_1 = [T_old_1(1);h_3_1+T_old_1(2);h_2_1+T_old_1(3);h_1_1+T_old_1(4);T_old_1(5);T_old_1(6);T_old_1(7);T_old_1(8)]; %c matrix with heat gen consideration on cone nodes
    T_new_1 = A_1\c_1; %finds new temp matrix given A/c matrices and old temp values
    
    J_old_1 = max(T_old_1); %finds max value from old temps
    J_new_1 = max(T_new_1); %finds max value from new temps
    
    ea_1 = abs((J_new_1-J_old_1)/J_new_1); %calculates error between old and new max temp, to be compared with stopping error at beginning of next loop
    
    T_old_1 = T_new_1; %redefines old temp matrix as current new temp matrix, recalculates new-new temps given old-new temps
    
    w_1(i) = J_new_1; %these three lines compile max temp data from each while loop run to plot it on one graph after stopping
    o_1(i) = p_1;
    i=i+1;
end

plot(o_1,w_1) %plots the compiled max temp data
title('Case 1 (q_s_h_e_a_r = 0 * q_f_r_i_c) Maximum Temperature Over Time')
xlabel('Time Step Number')
ylabel('Temperature [Celsius]')

E_in_1 = q_gen_1 * Area_cone * t_1;
E_out_1 = 0;

figure %creates new figure window for heatmap plot
Z_1 = [T_new_1(8);T_new_1(7);T_new_1(6);T_new_1(5);T_new_1(4);T_new_1(3);T_new_1(2)]; %flips order of T matrix to reflect having cutting tip on bottom in heatmap
H_1 = heatmap(Z_1, 'Title', 'Steady State Temperature Distribution, Case 1'); %creates a heatmap visual of the above Z matrix w/ color based on temp values
fprintf('-----------------------Case 1:--------------------------------\n\n') %readout title
fprintf('Time to SS: %.2f seconds\n', t_1) %prints time to steady state
fprintf('Highest Nodal Temp: %.2f Celsius\n', J_new_1) %prints highest temp value in celsius
fprintf('Total Energy Into Drill Bit: %.2f Joules\n', E_in_1)
fprintf('Total Energy out of Drill Bit: %.2f Joules\n\n', E_out_1)
fprintf('-----------------------Case 2:--------------------------------\n\n')

%--------------------------------------------------------------------------
%CASE 2 (q_shear = 0.5 * q_fric):

figure
while ea_2 > es
    p_2 = p_2 + 1; %time step counter
    t_2 = p_2 * del_t; %current time value
    c_2 = [T_old_2(1);h_3_2+T_old_2(2);h_2_2+T_old_2(3);h_1_2+T_old_2(4);T_old_2(5);T_old_2(6);T_old_2(7);T_old_2(8)];
    T_new_2 = A_2\c_2; %finds new temp matrix given A/c matrices and old temp values
    
    J_old_2 = max(T_old_2); %finds max value from old temps
    J_new_2 = max(T_new_2); %finds max value from new temps
    
    ea_2 = abs((J_new_2-J_old_2)/J_new_2); %calculates error between old and new max temp, to be compared with stopping error at beginning of next loop
    
    T_old_2 = T_new_2; %redefines old temp matrix as current new temp matrix, recalculates new-new temps given old-new temps
    
    w_2(j) = J_new_2;
    o_2(j) = p_2;
    j=j+1;
end

plot(o_2,w_2) 
title('Case 2 (q_s_h_e_a_r = 0.5 * q_f_r_i_c) Maximum Temperature Over Time')
xlabel('Time Step Number')
ylabel('Temperature [Celsius]')

E_in_2 = q_gen_2 * Area_cone * t_2;
E_out_2 = 0;

figure %makes a new figure window so Matlab doesn't replace last heatmap with this one (i.e. lets all three heatmaps display at once)
Z_2 = [T_new_2(8);T_new_2(7);T_new_2(6);T_new_2(5);T_new_2(4);T_new_2(3);T_new_2(2)];
H_2 = heatmap(Z_2, 'Title', 'Steady State Temperature Distribution, Case 2'); %creates a heatmap visual of the above Z matrix w/ color based on temp values
fprintf('Time to SS: %.2f seconds\n', t_2)
fprintf('Highest Nodal Temp: %.2f Celsius\n', J_new_2)
fprintf('Total Energy Into Drill Bit: %.2f Joules\n', E_in_2)
fprintf('Total Energy out of Drill Bit: %.2f Joules\n\n', E_out_2)
fprintf('-----------------------Case 3:--------------------------------\n\n')

%--------------------------------------------------------------------------
%CASE 3 (q_shear = 1 * q_fric):

figure
while ea_3 > es
    p_3 = p_3 + 1; %time step counter
    t_3 = p_3 * del_t; %current time value
    c_3 = [T_old_3(1);h_3_3+T_old_3(2);h_2_3+T_old_3(3);h_1_3+T_old_3(4);T_old_3(5);T_old_3(6);T_old_3(7);T_old_3(8)];
    T_new_3 = A_3\c_3; %finds new temp matrix given A/c matrices and old temp values
    
    J_old_3 = max(T_old_3); %finds max value from old temps
    J_new_3 = max(T_new_3); %finds max value from new temps
    
    ea_3 = abs((J_new_3-J_old_3)/J_new_3); %calculates error between old and new max temp, to be compared with stopping error at beginning of next loop
    
    T_old_3 = T_new_3; %redefines old temp matrix as current new temp matrix, recalculates new-new temps given old-new temps
    
    w_3(k) = J_new_3;
    o_3(k) = p_3;
    k=k+1;
end

plot(o_3,w_3)
title('Case 3 (q_s_h_e_a_r = 1 * q_f_r_i_c) Maximum Temperature Over Time')
xlabel('Time Step Number')
ylabel('Temperature [Celsius]')

E_in_3 = q_gen_3 * Area_cone * t_3;
E_out_3 = 0;

figure %makes a new figure window so Matlab doesn't replace last heatmap with this one (i.e. lets all three heatmaps display at once)
Z_3 = [T_new_3(8);T_new_3(7);T_new_3(6);T_new_3(5);T_new_3(5);T_new_3(3);T_new_3(2)];
H_3 = heatmap(Z_3, 'Title', 'Steady State Temperature Distribution, Case 3'); %creates a heatmap visual of the above Z matrix w/ color based on temp values
fprintf('Time to SS: %.2f seconds\n', t_3)
fprintf('Highest Nodal Temp: %.2f Celsius\n', J_new_3)
fprintf('Total Energy Into Drill Bit: %.2f Joules\n', E_in_3)
fprintf('Total Energy out of Drill Bit: %.2f Joules\n', E_out_3)