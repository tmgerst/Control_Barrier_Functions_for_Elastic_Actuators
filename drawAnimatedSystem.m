function [fig_simulation,video_handle] = drawAnimatedSystem(t,x,x_ref,sys_params)
%
    % Load system parameters
    l1 = sys_params.l1; l2 = sys_params.l2;
    l = l1+l2;
    
    % Compute joint positions upfront
    x_pos_link1 = l1*cos(x(:,1)); 
    y_pos_link1 = l1*sin(x(:,1)); 
    
    x_pos_link2 = x_pos_link1 + l2*cos(x(:,1)+x(:,5));
    y_pos_link2 = y_pos_link1 + l2*sin(x(:,1)+x(:,5));

    manipulator_positions.x = x_pos_link2;
    manipulator_positions.y = y_pos_link2;
    
    video_handle = VideoWriter('animation');
    video_handle.FrameRate = 15;
    open(video_handle);
    
    fig_simulation = figure('Name', 'Simulation');
    fig_simulation.Color = [1 1 1];
    labels = [];
    
    % Initialize objects
    txt1 = text(0.9,0.9,"");  txt1.Units = "normalized";
    txt2 = text(0.9,0.8,"");  txt2.Units = "normalized";
    txt3 = text(0.9,0.7,"");  txt3.Units = "normalized";

    % Animate system
    tic;
    i = 0; increment = 20;
    while i + increment < length(t)   
        i = i + increment;
        
        % Draw robot at current position
        figure(fig_simulation); clf(fig_simulation);

        scatter(0, 0, 50, 'c', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; % Base Point / First Joint
        scatter(x_pos_link1(i), y_pos_link1(i), 50, 'c', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Second Joint
        scatter(x_pos_link2(i), y_pos_link2(i), 50, 'r', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Manipulator

        xlim([-4/3*l, 4/3*l]); ylim([-4/3*l, 4/3*l]);
        
        plot([0, x_pos_link1(i)], [0, y_pos_link1(i)], 'Color', [0.75,0.75,0.75], 'LineWidth', 2, 'HandleVisibility', 'off'); % First Link
        plot([x_pos_link1(i), x_pos_link2(i)], [y_pos_link1(i), y_pos_link2(i)], 'Color', [0.75,0.75,0.75], 'LineWidth', 2, 'HandleVisibility', 'off'); % Second Link
        plot(manipulator_positions.x(1:i), manipulator_positions.y(1:i), 'g', 'LineWidth', 0.8);
        
        if i <= increment
            labels = [labels; "simulated trajectory"];
        end
        
        % Reference trajectory
        plot(l1*cos(x_ref.q1) + l2*cos(x_ref.q1+x_ref.q2), l1*sin(x_ref.q1) + l2*sin(x_ref.q1+x_ref.q2), ...
            'r--','LineWidth', 0.3);
        
        if i <= increment
            labels = [labels; "reference trajectory"];
        end
        
        % Constraint plotting
        x_pos_constraint_link1 = 1/3*l1*cos(3/4*pi); y_pos_constraint_link1 = 1/3*l1*sin(3/4*pi);
        
        x_pos_constraint_link2_upper = x_pos_link1(i) + 1/3*l2*cos(1/2*pi+x(i,1)); 
        x_pos_constraint_link2_lower = x_pos_link1(i) - 1/3*l2*cos(1/2*pi+x(i,1));
        y_pos_constraint_link2_upper = y_pos_link1(i) + 1/3*l2*sin(1/2*pi+x(i,1));
        y_pos_constraint_link2_lower = y_pos_link1(i) - 1/3*l2*sin(1/2*pi+x(i,1));
        
        plot([0, x_pos_constraint_link1], [0, y_pos_constraint_link1], 'k'); 
        plot([0, x_pos_constraint_link1], [0, -y_pos_constraint_link1], 'k');
        plot([x_pos_link1(i), x_pos_constraint_link2_upper], [y_pos_link1(i), y_pos_constraint_link2_upper], 'k'); 
        plot([x_pos_link1(i), x_pos_constraint_link2_lower], [y_pos_link1(i), y_pos_constraint_link2_lower], 'k');
        
        if i <= increment
            labels = [labels; "constraints"];
        end
        
        title("Animated motion of a 2-dof robot");
%         title('Planar robot with 2 degrees of freedom');

        xlabel("[m]"); ylabel("[m]");
        legend(labels,'Location','northwest');
        grid on;

        % Draw variable values
        txt1 = text(0.8, 0.95, "t="+t(i)+"s", 'Units', 'normalized');
        txt2 = text(0.8, 0.85, "q_1="+x(i,1)/pi*180+"°", 'Units', 'normalized'); hold off;
        txt3 = text(0.8, 0.80, "q_2="+x(i,5)/pi*180+"°", 'Units', 'normalized');
        
        frame = getframe(gcf);
        writeVideo(video_handle,frame);
    end
    toc;
    
    close(video_handle);
end