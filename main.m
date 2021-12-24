%
% 05.11.2021

clear all;
close all;
clc;

%% Parameters for the system
sys_params.l1 = 1;      sys_params.l2 = 1;
sys_params.m1 = 1;      sys_params.m2 = 1;
sys_params.r1 = 1;      sys_params.r2 = 1;
sys_params.k1 = 1e4;    sys_params.k2 = 1e4;
sys_params.d1 = 0;      sys_params.d2 = 0;
sys_params.j = 1;

draw_reference_states   = false;
draw_system_states      = false;
animate_system          = true;
draw_trajectory         = false;

sample_time = 0.01;
stop_time = 30;
t_span = (0:sample_time:stop_time)';

%% Construct reference trajectory
% xpos_ref = linspace(0.4,0,stop_time*1/sample_time)';
% % xpos_ref = [xpos_ref; -flip(xpos_ref)];
% ypos_ref = 0.5+0.1*sin(16*pi*xpos_ref);
% [theta1, theta2] = twoDofInverseKinematics(xpos_ref,ypos_ref,sys_params.l1,sys_params.l2);

x_ref = constructReferenceStates([],[],t_span);

bounds_min = [-3/4*pi;NaN;NaN;NaN;
              -pi/2;NaN;NaN;NaN];
bounds_max = [3/4*pi;NaN;NaN;NaN;
              pi/2;NaN;NaN;NaN];
[barrier_functions,grad_barrier_functions] = constructBarrierFunctions(bounds_min,bounds_max);
if isempty(barrier_functions)
    fprintf("No barrier functions have been constructed. Resuming simulation without constraint enforcement by CBFs.\n");
end

if draw_reference_states
    fig_reference_states = drawSystemStates(t_span,[],x_ref,bounds_min,bounds_max,'Reference');
end

%% System without damping
x0 = [x_ref.theta1(1); x_ref.theta1_first_dv(1); x_ref.theta1_second_dv(1); x_ref.theta1_third_dv(1);
      x_ref.theta2(1); x_ref.theta2_first_dv(1); x_ref.theta1_second_dv(1); x_ref.theta1_third_dv(1)];
% [t, x] = ode45(@(t,x) oneDofPlanarRobot(t,x,sys_params,x_ref),t_span,x0);
[t, x] = ode45(@(t,x) twoDofPlanarRobot(t,x,sys_params,x_ref,barrier_functions,grad_barrier_functions), t_span, x0);

 %% System with damping
% x0 = [x_ref.theta1(1), x_ref.theta1_first_dv(1), 0, 0, ...
%       x_ref.theta2(1), x_ref.theta2_first_dv(1), 0, 0];
% 
% var_names = {'t', 'q1', 'q1_dot', 'theta1', 'theta1_dot', 'q2', 'q2_dot', 'theta2', 'theta2_dot', ...
%              'v1', 'v2', 'u_dot1', 'u_dot2', 'u1', 'u2', 'tau1', 'tau2'};
% var_types = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
%             'double', 'double', 'double', 'double', 'double', 'double' 'double', 'double'}; 
% system_states = table('Size', [length(t_span),length(var_names)], 'VariableNames', var_names, 'VariableTypes', var_types);
% 
% for i=1:(length(t_span)-1)
%     %% Start conditions
%     if i==1
%         current_x = x0;
%         u = zeros(2,1);
%         tau = zeros(2,1);
%     else
%         current_x = x(end,:);
%     end
%     
%     %% Control via dynamic compensator
%     [v,index] = controllerForDamping(t_span,i,current_x,x_ref,sys_params);
%     u_dot = dynamicCompensator(current_x,v,u,sys_params);
%     
%     if i~=1
%         u = [trapz([system_states.t(1:i-1);t_span(i)],[system_states.u_dot1(1:i-1);u_dot(1)]); 
%             trapz([system_states.t(1:i-1);t_span(i)],[system_states.u_dot2(1:i-1);u_dot(2)])];
%     end
%     
%     tau = sys_params.j*eye(2)*u + ...
%         [sys_params.d1,0;0,sys_params.d2] * ([current_x(4);current_x(8)]-[current_x(2);current_x(6)]) + ...
%         [sys_params.k1,0;0,sys_params.k2] * ([current_x(3);current_x(7)]-[current_x(1);current_x(5)]);  
%     
%     %% Enforce constraints via CBFs
%     [q_second_deriv,~,~,theta_second_deriv] = StateVariablesHigherDerivatives(current_x,sys_params);
%     
%     B_max = pi/2 + 1/sys_params.k1 * ( tau(1)-sys_params.d1*(current_x(4)-current_x(2))-sys_params.j*q_second_deriv(1) ) - current_x(3);
%     grad_B_max = [0; sys_params.d1/sys_params.k1; -1; -sys_params.d1/sys_params.k1;...
%                   0; 0; 0; 0];
%     
%     f_x = [current_x(2);q_second_deriv(1);current_x(4);theta_second_deriv(1);...
%            current_x(6);q_second_deriv(2);current_x(8);theta_second_deriv(2)];
%               
%     inv_motor_inertia = inv(sys_params.j*eye(2));
%     g_x = [0, 0; 0, 0; 0, 0; inv_motor_inertia(1,1), inv_motor_inertia(1,2);
%            0, 0; 0, 0; 0, 0; inv_motor_inertia(2,1), inv_motor_inertia(2,2)];
%     
%     A = -[dot(grad_B_max,g_x(:,1)), dot(grad_B_max,g_x(:,2))];
%     b = [dot(grad_B_max,f_x)+B_max];
%     
%     opts_quadprog = optimset('Display','none');
%     tau = quadprog(eye(2),-tau',A,b,[],[],[],[],[],opts_quadprog);
%     
%     %% Save system states
%     system_states = copySystemStatesToTable(system_states,i,t_span(i),current_x,v,u_dot,u,tau);
%     
%     %% Compute system at next time step
%     [t,x] = ode45(@(t,x) twoDofPlanarRobotWithDamping(t,x,tau,sys_params), [t_span(i),t_span(i+1)], current_x);
%     fprintf("%i | %.4f\n", index, t(end));
% end

%% Plots and Animation
if draw_system_states
    fig_system_states = drawSystemStates(t,x,x_ref,bounds_min,bounds_max,'System');
end

if animate_system
    [fig_simulation, video_handle] = drawAnimatedSystem(t,x,x_ref,sys_params);
end

