function [theta1, theta2] = twoDofInverseKinematics(x,y,l1,l2)
%
    theta2 = acos((x.^2+y.^2-l1^2-l2^2)/(2*l1*l2));
    theta1 = atan(y./x) - atan(l2*sin(theta2)./(l1+l2*cos(theta2)));
    
    index = x<0;
    theta1(index) = theta1(index)+pi;
    theta1(3001) = (theta1(3000)+theta1(3002))/2;
    
    figure;
    plot(l1*(cos(theta1))+l2*cos(theta1+theta2), l1*sin(theta1)+l2*sin(theta1+theta2)); hold on;
    plot(x,y);
    legend('Inverse kinematics trajectory', 'Normal trajectory');
end