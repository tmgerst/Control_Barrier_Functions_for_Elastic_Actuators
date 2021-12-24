function [B, grad_B] = constructBarrierFunctions(bounds_min,bounds_max)
%
        gamma = 16;
        beta = 0;
        grad_B_min = [];    grad_B_max = [];        
        
        %% CBFs for min. bounds
        counter_min = 0;
        if sum(isnan(bounds_min)) ~= length(bounds_min)
            if ~isnan(bounds_min(1))
                counter_min = counter_min+1;
                B_min(counter_min,1) = {@(x) -gamma^3*bounds_min(1)+gamma^3*x(1)+3*gamma^2*x(2)+3*gamma*x(3)+x(4)-beta};
                grad_B_min = [grad_B_min, [gamma^3; 3*gamma^2; 3*gamma; 1; 0; 0; 0; 0]];
            end
            
            if ~isnan(bounds_min(2))
                counter_min = counter_min+1;
                B_min(counter_min,1) = {@(x) -gamma^2*bounds_min(2)+gamma^2*x(2)+2*gamma*x(3)+x(4)};
                grad_B_min = [grad_B_min, [0; gamma^2; 2*gamma; 1; 0; 0; 0; 0]];                
            end
            
            if ~isnan(bounds_min(5))
                counter_min = counter_min+1;
                B_min(counter_min,1) = {@(x) -gamma^3*bounds_min(5)+gamma^3*x(5)+3*gamma^2*x(6)+3*gamma*x(7)-x(8)-beta};
                grad_B_min = [grad_B_min, [0; 0; 0; 0; gamma^3; 3*gamma^2; 3*gamma; 1]];
            end
            
            if ~isnan(bounds_min(6))
                counter_min = counter_min+1;
                B_min(counter_min,1) = {@(x) -gamma^2*bounds_min(6)+gamma^2*x(6)+2*gamma*x(7)+x(8)};
                grad_B_min = [grad_B_min, [0; 0; 0; 0; 0; gamma^2; 2*gamma; 1]];                
            end            
        end
        
        %% CBFs for max. bounds
        counter_max = 0;
        if sum(isnan(bounds_max)) ~= length(bounds_max)
            if ~isnan(bounds_max(1))
                counter_max = counter_max+1;
                B_max(counter_max,1) = {@(x) gamma^3*bounds_max(1)-gamma^3*x(1)-3*gamma^2*x(2)-3*gamma*x(3)-x(4)-beta};
                grad_B_max = [grad_B_max, [-gamma^3; -3*gamma^2; -3*gamma; -1; 0; 0; 0; 0]];
            end
            
            if ~isnan(bounds_max(2))
                counter_max = counter_max+1;
                B_max(counter_max,1) = {@(x) gamma^2*bounds_max(2)-gamma^2*x(2)-2*gamma*x(3)-x(4)};
                grad_B_max = [grad_B_max, [0; -gamma^2; -2*gamma; -1; 0; 0; 0; 0]];                
            end
            
            if ~isnan(bounds_max(5))
                counter_max = counter_max+1;
                B_max(counter_max,1) = {@(x) gamma^3*bounds_max(5)-gamma^3*x(5)-3*gamma^2*x(6)-3*gamma*x(7)-x(8)-beta};
                grad_B_max = [grad_B_max, [0; 0; 0; 0; -gamma^3; -3*gamma^2; -3*gamma; -1]];
            end
            
            if ~isnan(bounds_max(6))
                counter_max = counter_max+1;
                B_max(counter_max,1) = {@(x) gamma^2*bounds_max(6)-gamma^2*x(6)-2*gamma*x(7)-x(8)};
                grad_B_max = [grad_B_max, [0; 0; 0; 0; 0; -gamma^2; -2*gamma; -1]];                
            end            
        end
        
        if counter_min & counter_max
            B = [B_min; B_max];
            grad_B = [grad_B_min, grad_B_max];
            return;
        elseif counter_max
            B = B_max;
            grad_B = grad_B_max;
            return
        elseif counter_min
            B = B_min;
            grad_B = grad_B_min;
            return
        else
            B = [];
            grad_B = [];
            fprintf("Warning: State constraints are missing.\n");
            return
        end            
end