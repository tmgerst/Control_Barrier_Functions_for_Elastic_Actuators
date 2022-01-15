function u_dot = dynamicCompensator(x,v,u,tau,tau_first_deriv,sys_params)
%
    %% System parameters and convenient variable naming
    m1 = sys_params.m1; m2 = sys_params.m2;
    l1 = sys_params.l1; l2 = sys_params.l2;
    r1 = sys_params.r1; r2 = sys_params.r2;
    k1 = sys_params.k1; k2 = sys_params.k2;
    d1 = sys_params.d1; d2 = sys_params.d2;
    j = sys_params.j;
    
    g = 9.81;
    lc1 = l1/2; lc2 = l2/2;

    I1 = 1/4*m1*r1^2 + 1/12*m1*l1^2;
    I2 = 1/4*m2*r2^2 + 1/12*m2*l2^2;    

    %% Calculate second derivative of link angles 
    D = [d1, 0; 0, d2];
    K = [k1, 0; 0, k2];
    
    % Inertia
    M = [m1*lc1^2 + I1 + m2*(l1^2+lc2^2+2*l1*lc2*cos(x(5))) + I2, m2*(lc2^2+l1*lc2*cos(x(5))) + I2; 
                m2*(lc2^2+l1*lc2*cos(x(5))) + I2,                       m2*lc2^2 + I2           ];     
    
    %% Calculate higher derivatives of state variables
    [~,~,~,~,delta,~,~,~,~,~] = StateVariablesHigherDerivatives(x,tau,tau_first_deriv,sys_params);
        
    %% Calculate u_dot
    u_dot = inv(D)*(M*v + delta - K*u);
end