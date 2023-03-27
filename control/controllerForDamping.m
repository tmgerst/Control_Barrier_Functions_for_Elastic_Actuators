function [v, index] = controllerForDamping(current_t,current_x,tau,tau_first_deriv,x_ref,sys_params)
%
    % Feedback gain matrix
    L = [1e4, 1e3, 300, 10, 0, 0, 0, 0;
         0, 0, 0, 0, 1e4, 1e3, 300, 10];
     
    % Find timestamp corresponding reference trajectory values
    index = find(x_ref.t==current_t);
    
    x_ref_link_angles = [x_ref.q1(index),x_ref.q1_first_deriv(index),x_ref.q1_second_deriv(index),x_ref.q1_third_deriv(index),...
                         x_ref.q2(index),x_ref.q2_first_deriv(index),x_ref.q2_second_deriv(index),x_ref.q2_third_deriv(index)];
    
    % Compute higher derivatives at current timestamp to compare them against the reference
    [q_second_deriv,q_third_deriv,~,~,~,~,~,~,~,~] = stateVariablesHigherDerivatives(current_x,tau,tau_first_deriv,sys_params); 
    
    x_link_angles = [current_x(1),current_x(2),q_second_deriv(1),q_third_deriv(1),...
                     current_x(5),current_x(6),q_second_deriv(2),q_third_deriv(2)];
    
    % Compute controller output
    v = [x_ref.q1_fourth_deriv(index); x_ref.q2_fourth_deriv(index)] + L*(x_ref_link_angles - x_link_angles)';
%     v = L*(x_ref_link_angles - x_link_angles)';
end