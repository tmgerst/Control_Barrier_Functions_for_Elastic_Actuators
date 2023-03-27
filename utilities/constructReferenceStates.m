function x_ref = constructReferenceStates(q1, q2, t_span)
%
    if isempty(q1) & isempty(q2)
        b1 = pi/15;
        b2 = pi/15;
        q1 = pi*sin(b1*t_span);
        q2 = pi*sin(b2*t_span);
    end

    % Compeltely horizontal position for testing 
%     q1            = zeros(size(t_span));
%     q1_first_deriv   = zeros(size(t_span));
%     q1_second_deriv  = zeros(size(t_span));
%     q1_third_deriv   = zeros(size(t_span));
%     q1_fourth_deriv  = zeros(size(t_span));
%     
%     q2            = zeros(size(t_span));
%     q2_first_deriv   = zeros(size(t_span));
%     q2_second_deriv  = zeros(size(t_span));
%     q2_third_deriv   = zeros(size(t_span));
%     q2_fourth_deriv  = zeros(size(t_span));
    
    q1_first_deriv     = pi*b1*cos(b1*t_span);
    q1_second_deriv    = -pi*b1^2*sin(b1*t_span);
    q1_third_deriv     = -pi*b1^3*cos(b1*t_span);
    q1_fourth_deriv    = pi*b1^4*sin(b1*t_span);
      
    q2_first_deriv     = pi*b2*cos(b2*t_span);
    q2_second_deriv    = -pi*b2^2*sin(b2*t_span);
    q2_third_deriv     = -pi*b2^3*cos(b2*t_span);
    q2_fourth_deriv    = pi*b2^4*sin(b2*t_span);
    
    t = t_span;
    x_ref = table(t, ...
        q1, q1_first_deriv, q1_second_deriv, q1_third_deriv, q1_fourth_deriv, ...
        q2, q2_first_deriv, q2_second_deriv, q2_third_deriv, q2_fourth_deriv);
end