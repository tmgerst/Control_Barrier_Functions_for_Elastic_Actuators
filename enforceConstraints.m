function u_enforced = enforceConstraints(control_for_damping_used,current_x,sys_params,u_nom,tau,tau_first_deriv,barrier_functions,grad_barrier_functions)
%       
        if isempty(barrier_functions)
            u_enforced = u_nom;
            return;
        end
               
        %% Enforcement through inequalities (experimental, does not work well)
%         f_x_qp_help =  repmat(f_x_qp,1,size(grad_barrier_functions,2));
%         g_x_qp_help1 = repmat(g_x_qp1,1,size(grad_barrier_functions,2));
%         g_x_qp_help2 = repmat(g_x_qp2,1,size(grad_barrier_functions,2));        
%         
%         A = [dot(grad_barrier_functions,g_x_qp_help1)', dot(grad_barrier_functions,g_x_qp_help2)'];
%         b = [dot(grad_barrier_functions,f_x_qp_help)'+cellfun(@(v) v(x),barrier_functions)];     
%         
% %         lower_constraints_admissable = A(1,:)*u_nom+b(1) >= 0;
% %         upper_constraints_admissable = A(2,:)*u_nom+b(2) >= 0;
%         
%         lower_constraints_admissable = A(1:2,:)*u_nom+b(1:2) >= 0;
%         upper_constraints_admissable = A(3:4,:)*u_nom+b(3:4) >= 0;
%         
%         u_enforced = zeros(2,1);
%         if all(lower_constraints_admissable) & all(upper_constraints_admissable)
%             u_enforced = u_nom;
%         end
%         if ~all(lower_constraints_admissable)
% %             u_enforced = A(1,:)\(-b(1));
%             u_enforced = A(1:2,:)\(-b(1:2));
% 
% %             if ~lower_constraints_admissable(1)
% %                 u_enforced = [u_enforced(1);u_nom(2)];
% %             else
% %                 u_enforced = [u_nom(1);u_enforced(2)];
% %             end
%         end
%         if ~all(upper_constraints_admissable)
% %             u_enforced = A(2,:)\(-b(2));
%             u_enforced = A(3:4,:)\(-b(3:4));
%             
% %             if ~upper_constraints_admissable(1)
% %                 u_enforced = [u_enforced(1);u_nom(2)];
% %             else
% %                 u_enforced = [u_nom(1);u_enforced(2)];
% %             end
%         end  
        
        %% Quadratic program for damped model
        [q_second_deriv,q_third_deriv,~,~,~,M,M_dot,M_dotdot,n,n_dotdot] = ...
            StateVariablesHigherDerivatives(current_x,tau,tau_first_deriv,sys_params);
        
        if control_for_damping_used
            f_x = [current_x(2); q_second_deriv(1); q_third_deriv(1); 0;
                   current_x(6); q_second_deriv(2); q_third_deriv(2); 0];

            g_x = [0, 0; 0, 0; 0, 0; 1, 0;
                   0, 0; 0, 0; 0, 0; 0, 1];
        else
            K = [sys_params.k1, 0; 0, sys_params.k2];
            J = sys_params.j * eye(2);
            
            f_x4 = -inv(M)*K*inv(J)*(M*q_second_deriv+n)...
                -inv(M)*( M_dotdot*q_second_deriv+2*M_dot*q_third_deriv+n_dotdot+K*q_second_deriv );

            g_x4 = inv(M) * K * inv(J);
            
            f_x = [current_x(2); q_second_deriv(1); q_third_deriv(1); f_x4(1);
                   current_x(6); q_second_deriv(2); q_third_deriv(2); f_x4(2)];
            g_x = [0, 0; 0, 0; 0, 0; g_x4(1,1), g_x4(1,2);
                   0, 0; 0, 0; 0, 0; g_x4(2,1), g_x4(2,2)];
        end
        
        f_x_qp =  repmat(f_x,1,size(grad_barrier_functions,2));
        g_x_qp1 = repmat(g_x(:,1),1,size(grad_barrier_functions,2));
        g_x_qp2 = repmat(g_x(:,2),1,size(grad_barrier_functions,2));
        
        A = -[dot(grad_barrier_functions,g_x_qp1)', dot(grad_barrier_functions,g_x_qp2)'];
        
        % Barrier function input has structure: 
        % [q1, q1_first_deriv, q1_second_deriv, q1_third_deriv, q2, q2_first_deriv, q2_second_deriv, q2_third_deriv]
        barrier_function_input = [current_x(1), current_x(2), q_second_deriv(1), q_third_deriv(1),...
                                  current_x(5), current_x(6), q_second_deriv(2), q_third_deriv(2)];
        b = [dot(grad_barrier_functions,f_x_qp)'+cellfun(@(v) v(barrier_function_input),barrier_functions)];
        
        Aeq = [];   beq = [];
        lb = [];    ub = [];
        nonlcon = [];
        
        opts_quadprog = optimset('Display','none');
%         opts_quadprog.StepTolerance = 1e-2; 
%         opts_quadprog.OptimalityTolerance = 1e-2;
        u_enforced = quadprog(eye(2),-u_nom',A,b,[],[],[],[],[],opts_quadprog);
        
        if isempty(u_enforced)
            dbstop at 80;
        end
end