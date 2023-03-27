function u_nom = systemController(t,x,x_ref,f_x,g_x,dof_degree)
%% IMPORTANT: Legacy function, needs to be taken out of the repo
%
    if dof_degree == 1
        %% 1-dof-system: Control Input
        theta1              = interp1(x_ref.t,x_ref.theta1,t);
        theta1_first_dv     = interp1(x_ref.t,x_ref.theta1_first_dv,t);
        theta1_second_dv    = interp1(x_ref.t,x_ref.theta1_second_dv,t);
        theta1_third_dv     = interp1(x_ref.t,x_ref.theta1_third_dv,t);
        theta1_fourth_dv    = interp1(x_ref.t,x_ref.theta1_fourth_dv,t);
        
        L = [30, 20, 10, 2];
        v = theta1_fourth_dv + L*([theta1; theta1_first_dv; theta1_second_dv; theta1_third_dv] - x);
        u_nom = inv(g_x)*(v-f_x);
        
        %% 1-dof system: CBFs and quadratic program
        gamma = 1; % tuning parameter, an increase leads to the system getting closer to its boundary
        boundary_theta1_max = 1/2*pi;
        boundary_theta1_min = -1/2*pi;
        
        B_x_theta1_max = @(x) gamma^3*boundary_theta1_max-gamma^3*x(1)-3*gamma^2*x(2)-3*gamma*x(3)-x(4);
        grad_B_theta1_max = [-gamma^3; -3*gamma^2; -3*gamma; -1];
        
        B_x_theta1_min = @(x) -gamma^3*boundary_theta1_min+gamma^3*x(1)+3*gamma^2*x(2)+3*gamma*x(3)+x(4);
        grad_B_theta1_min = [gamma^3; 3*gamma^2; 3*gamma; 1];
        
        f_x_qp = [x(2); x(3); x(4); f_x];
        g_x_qp = [0; 0; 0; g_x];

        % Inequality constraints with -A(x)*u <= -b
        A = [dot(grad_B_theta1_max,g_x_qp);
             dot(grad_B_theta1_min,g_x_qp)];
        b = -[dot(grad_B_theta1_max,f_x_qp)+B_x_theta1_max(x); 
            dot(grad_B_theta1_min,f_x_qp)+B_x_theta1_min(x)];

        indices = b./A;
        if u_nom > indices(1)
            u_nom = indices(1);
        end
        if u_nom < indices(2)
            u_nom = indices(2);
        end
        
        fprintf("t: "+t+" | "+indices(1)+" | "+indices(2)+"\n");
        
%         fprintf("t: %.6f | Theta1: %.4f | u_nom-u: %.2f\n", t, x(1), u_nom(1)-u_nom(1));
        
    elseif dof_degree == 2
        %% 2-dof-system: Control Input
        theta1              = interp1(x_ref.t,x_ref.theta1,t);
        theta1_first_dv     = interp1(x_ref.t,x_ref.theta1_first_dv,t);
        theta1_second_dv    = interp1(x_ref.t,x_ref.theta1_second_dv,t);
        theta1_third_dv     = interp1(x_ref.t,x_ref.theta1_third_dv,t);
        theta1_fourth_dv    = interp1(x_ref.t,x_ref.theta1_fourth_dv,t);

        theta2              = interp1(x_ref.t,x_ref.theta2,t);
        theta2_first_dv     = interp1(x_ref.t,x_ref.theta2_first_dv,t);
        theta2_second_dv    = interp1(x_ref.t,x_ref.theta2_second_dv,t);
        theta2_third_dv     = interp1(x_ref.t,x_ref.theta2_third_dv,t);
        theta2_fourth_dv    = interp1(x_ref.t,x_ref.theta2_fourth_dv,t);

        L = [1e4, 1e3, 300, 10, 0, 0, 0, 0;
             0, 0, 0, 0, 1e4, 1e3, 300, 10];

        v =  [theta1_fourth_dv; theta2_fourth_dv] + L*( [theta1; theta1_first_dv; theta1_second_dv; theta1_third_dv;
                                                         theta2; theta2_first_dv; theta2_second_dv; theta2_third_dv] - x);

        u_nom = inv(g_x)*(v-f_x); 
    end
end