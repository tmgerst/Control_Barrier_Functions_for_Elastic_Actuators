function u_enforced = enforceConstraints(x,u_nom,f_x,g_x,barrier_functions,grad_barrier_functions)
%       
        if isempty(barrier_functions)
            u_enforced = u_nom;
            return;
        end

        f_x_qp = [x(2); x(3); x(4); f_x(1); 
                  x(6); x(7); x(8); f_x(2)];
        g_x_qp1 = [0; 0; 0; g_x(1,1);
                   0; 0; 0; g_x(2,1)];
        g_x_qp2 = [0; 0; 0; g_x(1,2);
                   0; 0; 0; g_x(2,2)];
               
        %% Enforcement through inequalities
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
        
        %% Quadratic program
        % Objective Function and start point
%         quadratic_program = @(u)dot(u,u)-2*dot(u_nom,u);
%         u0 = u_nom;
        
        f_x_qp_help =  repmat(f_x_qp,1,size(grad_barrier_functions,2));
        g_x_qp_help1 = repmat(g_x_qp1,1,size(grad_barrier_functions,2));
        g_x_qp_help2 = repmat(g_x_qp2,1,size(grad_barrier_functions,2));
        
        A = -[dot(grad_barrier_functions,g_x_qp_help1)', dot(grad_barrier_functions,g_x_qp_help2)'];
        b = [dot(grad_barrier_functions,f_x_qp_help)'+cellfun(@(v) v(x),barrier_functions)];
        
        Aeq = [];   beq = [];
        lb = [];    ub = [];
        nonlcon = [];

%         opts = optimoptions('fmincon','MaxFunctionEvaluations',100,'MaxIterations',50,...
%             'OptimalityTolerance',1e-2,'StepTolerance',1e-2,'FiniteDifferenceStepSize',1e-2,...
%             'Display','none');
%         u_enforced = fmincon(quadratic_program,u0,A,b,Aeq,beq,lb,ub,nonlcon,opts);
        
        opts_quadprog = optimset('Display','none');
        opts_quadprog.StepTolerance = 1e-2; 
        opts_quadprog.OptimalityTolerance = 1e-2;
        u_enforced = quadprog(eye(2),-u_nom',A,b,[],[],[],[],[],opts_quadprog);
end