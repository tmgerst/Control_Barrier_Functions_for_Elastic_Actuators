function [v, index] = controllerForDamping(t_span,i,current_x,x_ref,sys_params)
%
    L = [1e4, 1e3, 300, 10, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 1e4, 1e3, 300, 10, 0];

%     L = [1, 10, 0, 0, 0, 0, 0, 0, 0, 0;
%          0, 0, 0, 0, 0, 1, 10, 0, 0, 0];

    [q_second_deriv,q_third_deriv,q_fourth_deriv,~] = StateVariablesHigherDerivatives(current_x,sys_params); 
     
    index = find(x_ref.t==t_span(i));
    v = L*( [x_ref.theta1(index),x_ref.theta1_first_dv(index),x_ref.theta1_second_dv(index),x_ref.theta1_third_dv(index),x_ref.theta1_fourth_dv(index),...
             x_ref.theta2(index),x_ref.theta2_first_dv(index),x_ref.theta2_second_dv(index),x_ref.theta2_third_dv(index),x_ref.theta2_fourth_dv(index)] - ...
            [current_x(1),current_x(2),q_second_deriv(1),q_third_deriv(1),q_fourth_deriv(1),...
             current_x(5),current_x(6),q_second_deriv(2),q_third_deriv(2),q_fourth_deriv(2)] )';
end