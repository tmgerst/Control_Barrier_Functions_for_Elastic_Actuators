function fig_states = drawSystemStates(t,x,x_ref,bounds_min,bounds_max,decision_flag)
%
    if decision_flag == "System"
        fig_states = figure('Name', 'System states');
        fig_states.Color = [1 1 1];
        titles = {'$q_1 [rad]$', '$q_2 [rad]$', '$\dot{q_1} [rad/s]$', '$\dot{q_2} [rad/s]$', ...
                 '$\ddot{q_1} \ [rad/s^2]$', '$\ddot{q_2} \ [rad/s^2]$', '$q_1^{(3)} \ [rad/s^3]$', '$q_2^{(3)} \ [rad/s^3]$'};

        labels = [];
        counter = 0;
        for i=[1,3,5,7] % States for theta 1
            counter = counter + 1; 
            subplot(4,2,i); 
            plot(t, x(:,counter),'g','LineWidth',1.5); hold on; 
            plot(t, x_ref{:,counter+1},'r--');
            
            if i==1
                labels = [labels; "simulated trajectory"; "reference trajectory"];
            end
            
            if i==1
                if ~isnan(bounds_min(1))
                    plot(t,bounds_min(1)*ones(size(t)),'k','LineWidth',1.1);
                end
                if ~isnan(bounds_max(1))
                    plot(t,bounds_max(1)*ones(size(t)),'k','LineWidth',1.1);
                end
            end
            
            if i==3 
                if ~isnan(bounds_min(2))
                    plot(t,bounds_min(2)*ones(size(t)),'k','LineWidth',1.1);
                end
                if ~isnan(bounds_max(2))
                    plot(t,bounds_max(2)*ones(size(t)),'k','LineWidth',1.1);
                end
            end
            title(titles{i}, 'Interpreter', 'latex');
            xlabel('Time [s]'); grid on;
        end
        
        for i=[2,4,6,8]
            counter = counter + 1;
            subplot(4,2,i); 
            plot(t, x(:,counter),'g','LineWidth',1.5); hold on; plot(t, x_ref{:,counter+2},'r--');
            
            if i==2
                if ~isnan(bounds_min(5))
                    plot(t,bounds_min(5)*ones(size(t)),'k','LineWidth',1.1);
                    labels = [labels; "constraints"];
                end
                if ~isnan(bounds_max(5))
                    plot(t,bounds_max(5)*ones(size(t)),'k','LineWidth',1.1);
                end
            end
            
            if i==4
                if ~isnan(bounds_min(6))
                    plot(t,bounds_min(6)*ones(size(t)),'k','LineWidth',1.1);
                end
                if ~isnan(bounds_max(6))
                    plot(t,bounds_max(6)*ones(size(t)),'k','LineWidth',1.1);
                end
            end
            
            title(titles{i}, 'Interpreter', 'latex');
            xlabel('Time [s]'); grid on;
        end
        
        subplot(422);
        lgd = legend(labels); lgd.Position = [0.462760416666667,0.676729201573669,0.109114583333333,0.074094401756311];
        sgtitle('Link angles with derivatives');
        
    elseif decision_flag == "Reference"
        % x_ref ist eine Tabelle!
        fig_states = figure('Name', 'Reference states');
        titles = {'$\theta_1 [rad]$', '$\theta_2 [rad]$', '$\dot{\theta_1} [rad/s]$', '$\dot{\theta_2} [rad/s]$', ...
                 '$\ddot{\theta_1} \ [rad/s^2]$', '$\ddot{\theta_2} \ [rad/s^2]$', '$\theta_1^{(3)}$', '$\theta_2^{(3)}$', ...
                 '$\theta_1^{(4)}$', '$\theta_2^{(4)}$'};
        
        counter = 1;
        for i=[1,3,5,7,9]
            counter = counter + 1;
            subplot(5,2,i); plot(x_ref.t, x_ref{:,counter}); title(titles{i}, 'Interpreter', 'latex');
            xlabel('Time [s]'); grid on;
        end
        
        for i=[2,4,6,8,10]
            counter = counter + 1;
            subplot(5,2,i); plot(x_ref.t, x_ref{:,counter}); title(titles{i}, 'Interpreter', 'latex');
            xlabel('Time [s]'); grid on;
        end
        
        sgtitle('Reference states');
    end
end