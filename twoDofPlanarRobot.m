function dxdt = twoDofPlanarRobot(t,x,sys_params,x_ref,barrier_functions,grad_barrier_functions)
% x = [theta1, theta1_d, theta1_dd, theta1_ddd;
%      theta2, theta2_d, theta2_dd, theta2_ddd];

    %% Hilfsbenennungen
    dof_degree = 2;
    
    m1 = sys_params.m1; m2 = sys_params.m2;
    l1 = sys_params.l1; l2 = sys_params.l2;
    r1 = sys_params.r1; r2 = sys_params.r2;
    k1 = sys_params.k1; k2 = sys_params.k2;
    j = sys_params.j;
    
    g = 9.81;
    lc1 = l1/2; lc2 = l2/2;

    I1 = 1/4*m1*r1^2 + 1/12*m1*l1^2;
    I2 = 1/4*m2*r2^2 + 1/12*m2*l2^2;     
    
%     x(1) = wrapToPi(x(1));
%     x(5) = wrapToPi(x(5));
    
    %% Trägheitsmatrizen mit Ableitungen
    M = [m1*lc1^2 + I1 + m2*(l1^2+lc2^2+2*l1*lc2*cos(x(5))) + I2, m2*(lc2^2+l1*lc2*cos(x(5))) + I2; 
                    m2*(lc2^2+l1*lc2*cos(x(5))) + I2,                       m2*lc2^2 + I2           ];
                
    M_dot = [-2*m2*l1*lc2*sin(x(5))*x(6), -m2*l1*lc2*sin(x(5))*x(6);
             -m2*l1*lc2*sin(x(5))*x(6),                   0                   ];
    
    M_dotdot = [-2*m2*l1*lc2*(cos(x(5))*x(6)+sin(x(5))*x(7)), -m2*l1*lc2*(cos(x(5))*x(6)+sin(x(5))*x(7));
                -m2*l1*lc2*(cos(x(5))*x(6)+sin(x(5))*x(7)),                    0                        ];
    
    %% Coriolis-, Zentrifugalkräfte mit Ableitungen
    C = [(-m2*l1*lc2*sin(x(5)))*x(6)^2 - 2*m2*l1*lc2*sin(x(5))*x(2)*x(6); 
                            m2*l1*lc2*sin(x(5))*x(2)^2];
    
    c_expr1 = sin(x(5))*x(6)^4 - 3*cos(x(5))*x(6)^2*x(7);
    c_expr2 = -2*cos(x(5))*x(6)^2*x(7) - 2*sin(x(5))*(x(7)^2+x(6)*x(8));
    c_expr3 = 2*sin(x(5))*x(6)*x(7)*x(2) - 2*cos(x(5))*(x(3)*x(7)+x(2)*x(8));
    c_expr4 = -2*cos(x(5))*x(6)*(x(3)*x(6)+x(2)*x(7)) - 2*sin(x(5))*(x(4)*x(6)+2*x(3)*x(7)+x(2)*x(8));
    c_expr5 = -sin(x(5))*x(6)^2*x(2)^2 + cos(x(5))*(x(6)*x(2)^2+2*x(6)*x(2)*x(3));
    c_expr6 = 2*cos(x(5))*x(6)*x(2)*x(3) + 2*sin(x(5))*(x(3)^2+x(2)*x(4));
    
    C_dotdot = m2*l1*lc2 * [c_expr1 + c_expr2 + c_expr3 + c_expr4; 
                                c_expr5 + c_expr6];
    
    %% Gravitation mit Ableitungen
    G = [m1*lc1*g*cos(x(1)) + m2*g*(lc2*cos(x(1)+x(5))+l1*cos(x(1)));
                            m2*g*lc2*cos(x(1)+x(5))];
                        
    g_expr1 = -m1*lc2*g*cos(x(1))*x(2)^2 - m1*lc1*g*sin(x(1))*x(3);
    g_expr2 = -m2*g*(lc2*cos(x(1)+x(5))*(x(2)+x(6))^2+lc2*sin(x(1)+x(5))*(x(3)+x(7)));
    g_expr3 = -m2*g*(l1*cos(x(1))*x(2)^2+l1*sin(x(1))*x(3));
    g_expr4 = -m2*g*lc2*cos(x(1)+x(5))*(x(2)+x(6))^2;
    g_expr5 = -m2*g*lc2*sin(x(1)+x(5))*(x(3)+x(7));
    
    G_dotdot = [g_expr1 + g_expr2 + g_expr3;
                    g_expr4 + g_expr5];
    
    % C und G kombinieren analog Paper
    n = C + G;
    n_dotdot = C_dotdot + G_dotdot;
    
    %% Federmatrix und Trägheiten              
    K = [k1, 0;
         0, k2];
    
    J = j * eye(2);
    
    %% Einkoppelfunktion und Kontrollinput berechnen
    f_x = -inv(M)*K*inv(J)*(M*[x(3);x(7)]+n) - inv(M)*( M_dotdot*[x(3);x(7)]+2*M_dot*[x(4);x(8)]+n_dotdot+K*[x(3);x(7)] );
    g_x = inv(M)*K*inv(J);

    u_nom = systemController(t,x,x_ref,f_x,g_x,dof_degree);
    u = enforceConstraints(x,u_nom,f_x,g_x,barrier_functions,grad_barrier_functions);
    
    fprintf("t: %.6f | Theta1: %.4f | u_nom1-u1: %.2f | Theta2: %.4f | u_nom2-u2: %.2f\n", t, x(1), u_nom(1)-u(1), x(5), u_nom(2)-u(2));
    
    %% Nächsten Zustand berechnen
    % Vierte Ableitungen der jeweiligen Link-Winkel berechnen
    fourth_derivatives = f_x+g_x*u;
    
    dxdt(1) = x(2);
    dxdt(2) = x(3);
    dxdt(3) = x(4);
    dxdt(4) = fourth_derivatives(1);
    dxdt(5) = x(6);
    dxdt(6) = x(7);
    dxdt(7) = x(8);
    dxdt(8) = fourth_derivatives(2); 
    
    dxdt = dxdt';
end