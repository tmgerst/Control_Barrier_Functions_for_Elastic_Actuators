function x_ref = constructReferenceStates(theta1, theta2, t_span)
%
    if isempty(theta1) & isempty(theta2)
        b1 = pi/15;
        b2 = pi/15;
        theta1 = pi*sin(b1*t_span);
        theta2 = pi*sin(b2*t_span);
    end

    % Compeltely horizontal position for testing 
%     theta1            = zeros(size(t_span));
%     theta1_first_dv   = zeros(size(t_span));
%     theta1_second_dv  = zeros(size(t_span));
%     theta1_third_dv   = zeros(size(t_span));
%     theta1_fourth_dv  = zeros(size(t_span));
%     
%     theta2            = zeros(size(t_span));
%     theta2_first_dv   = zeros(size(t_span));
%     theta2_second_dv  = zeros(size(t_span));
%     theta2_third_dv   = zeros(size(t_span));
%     theta2_fourth_dv  = zeros(size(t_span));
    
    theta1_first_dv     = pi*b1*cos(b1*t_span);
    theta1_second_dv    = -pi*b1^2*sin(b1*t_span);
    theta1_third_dv     = -pi*b1^3*cos(b1*t_span);
    theta1_fourth_dv    = pi*b1^4*sin(b1*t_span);
      
    theta2_first_dv     = pi*b2*cos(b2*t_span);
    theta2_second_dv    = -pi*b2^2*sin(b2*t_span);
    theta2_third_dv     = -pi*b2^3*cos(b2*t_span);
    theta2_fourth_dv    = pi*b2^4*sin(b2*t_span);
    
    t = t_span;
    x_ref = table(t, ...
        theta1, theta1_first_dv, theta1_second_dv, theta1_third_dv, theta1_fourth_dv, ...
        theta2, theta2_first_dv, theta2_second_dv, theta2_third_dv, theta2_fourth_dv);
end