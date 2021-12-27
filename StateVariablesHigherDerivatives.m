function [q_second_deriv, q_third_deriv, q_fourth_deriv, theta_second_deriv, delta] = StateVariablesHigherDerivatives(x,tau,sys_params)
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
    J = j*eye(2);
    
    % Inertia
    M = [m1*lc1^2 + I1 + m2*(l1^2+lc2^2+2*l1*lc2*cos(x(5))) + I2, m2*(lc2^2+l1*lc2*cos(x(5))) + I2; 
                m2*(lc2^2+l1*lc2*cos(x(5))) + I2,                       m2*lc2^2 + I2           ];
    
    % Gravitational forces
    G = [m1*lc1*g*cos(x(1)) + m2*g*(lc2*cos(x(1)+x(5))+l1*cos(x(1)));
                            m2*g*lc2*cos(x(1)+x(5))];  
    
    % Coriolis forces
    C = [(-m2*l1*lc2*sin(x(5)))*x(6)^2 - 2*m2*l1*lc2*sin(x(5))*x(2)*x(6); 
                        m2*l1*lc2*sin(x(5))*x(2)^2];                    
    
    n = C + G;                              
    q_second_deriv = -inv(M) * ( n + D*([x(2);x(6)]-[x(4);x(8)]) + K*([x(1);x(5)]-[x(3);x(7)]) );

    %% Calculate third derivative of link angles
    M_dot = [-2*m2*l1*lc2*sin(x(5))*x(6), -m2*l1*lc2*sin(x(5))*x(6);
              -m2*l1*lc2*sin(x(5))*x(6),                   0                   ];
    
    % Coriolis forces first derivatives
    c_dot_expr1 = -cos(x(5))*x(5)^3;
    c_dot_expr2 = -2*sin(x(5))*x(6)*q_second_deriv(2);
    c_dot_expr3 = -2*cos(x(5))*x(2)*q_second_deriv(2);
    c_dot_expr4 = -2*sin(x(5))*(q_second_deriv(1)*x(6)+x(2)*q_second_deriv(2));
    c_dot_expr5 = cos(x(5))*x(6)*x(2)^2;
    c_dot_expr6 = 2*sin(x(5))*x(2)*q_second_deriv(1);
    
    C_dot = m2*l1*lc2*[c_dot_expr1+c_dot_expr2+c_dot_expr3+c_dot_expr4;
                                c_dot_expr5+c_dot_expr6];
    
    % Gravitational forces first derivatives
    g_dot_expr1 = -m1*lc1*g*sin(x(1))*x(2);
    g_dot_expr2 = -m2*g*lc2*sin(x(1)+x(5))*(x(2)+x(6));
    g_dot_expr3 = -m2*g*l1*sin(x(1))*x(2);
    g_dot_expr4 = -m2*g*lc2*sin(x(1)+x(5))*(x(2)+x(6));
    
    G_dot = [g_dot_expr1+g_dot_expr2+g_dot_expr3;
                    g_dot_expr4];
    
    n_dot = C_dot + G_dot;
    beta = n_dot + D*q_second_deriv + K*[x(2);x(6)] + M_dot*q_second_deriv;
    theta_second_deriv = inv(J)*tau - inv(J)*( D*([x(4);x(8)]-[x(2);x(6)]) + K*([x(3);x(7)]-[x(1);x(5)]) );
    
    q_third_deriv = inv(M) * ( -beta + D*theta_second_deriv + K*[x(4);x(8)] );
    
    %% Calculate fourth derivative of link angles
    M_dotdot = [-2*m2*l1*lc2*(cos(x(5))*x(6)+sin(x(5))*q_second_deriv(2)), -m2*l1*lc2*(cos(x(5))*x(6)+sin(x(5))*q_second_deriv(2));
                -m2*l1*lc2*(cos(x(5))*x(6)+sin(x(5))*q_second_deriv(2)),                    0                        ];
          
    % Coriolis forces second derivatives
    c_dotdot_expr1 = sin(x(5))*x(6)^4 - 3*cos(x(5))*x(6)^2*q_second_deriv(2);
    c_dotdot_expr2 = -2*cos(x(5))*x(6)^2*q_second_deriv(2) - 2*sin(x(5))*(q_second_deriv(2)^2+x(6)*q_third_deriv(2));
    c_dotdot_expr3 = 2*sin(x(5))*x(6)*q_second_deriv(2)*x(2) - 2*cos(x(5))*(q_second_deriv(1)*q_second_deriv(2)+x(2)*q_third_deriv(2));
    c_dotdot_expr4 = -2*cos(x(5))*x(6)*(q_second_deriv(1)*x(6)+x(2)*q_second_deriv(2)) - 2*sin(x(5))*(q_third_deriv(1)*x(6)+2*q_second_deriv(1)*q_second_deriv(2)+x(2)*q_third_deriv(2));
    c_dotdot_expr5 = -sin(x(5))*x(6)^2*x(2)^2 + cos(x(5))*(x(6)*x(2)^2+2*x(6)*x(2)*q_second_deriv(1));
    c_dotdot_expr6 = 2*cos(x(5))*x(6)*x(2)*q_second_deriv(1) + 2*sin(x(5))*(q_second_deriv(1)^2+x(2)*q_third_deriv(1));
    
    C_dotdot = m2*l1*lc2 * [c_dotdot_expr1 + c_dotdot_expr2 + c_dotdot_expr3 + c_dotdot_expr4; 
                                c_dotdot_expr5 + c_dotdot_expr6];
                        
    % Gravitational forces second derivatives
    g_dotdot_expr1 = -m1*lc2*g*cos(x(1))*x(2)^2 - m1*lc1*g*sin(x(1))*q_second_deriv(1);
    g_dotdot_expr2 = -m2*g*(lc2*cos(x(1)+x(5))*(x(2)+x(6))^2+lc2*sin(x(1)+x(5))*(q_second_deriv(1)+q_second_deriv(2)));
    g_dotdot_expr3 = -m2*g*(l1*cos(x(1))*x(2)^2+l1*sin(x(1))*q_second_deriv(1));
    g_dotdot_expr4 = -m2*g*lc2*cos(x(1)+x(5))*(x(2)+x(6))^2;
    g_dotdot_expr5 = -m2*g*lc2*sin(x(1)+x(5))*(q_second_deriv(1)+q_second_deriv(2));
    
    G_dotdot = [g_dotdot_expr1 + g_dotdot_expr2 + g_dotdot_expr3;
                    g_dotdot_expr4 + g_dotdot_expr5];
                
    n_dotdot = C_dotdot + G_dotdot;
       
    delta = n_dotdot + D*q_third_deriv + K*q_second_deriv + M_dotdot*q_second_deriv + 2*M_dot*q_third_deriv;   
    theta_third_deriv = -inv(J) * ( D*(theta_second_deriv-q_second_deriv) + K*([x(4);x(8)]-[x(2);x(6)]) );
    q_fourth_deriv = inv(M) * ( -delta + D*theta_third_deriv + K*theta_second_deriv);
end