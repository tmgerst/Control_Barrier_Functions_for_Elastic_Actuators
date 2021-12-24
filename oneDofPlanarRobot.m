function dxdt = oneDofPlanarRobot(t,x,sys_params,x_ref)
%
   dof_degree = 1;
    
   m1 = sys_params.m1;
   l1 = sys_params.l1;
   r1 = sys_params.r1;
   k1 = sys_params.k1;
   J = sys_params.j;
   
   g = 9.81;
   
   M = m1*(l1/2)^2 + 1/4*m1*r1^2+1/12*m1*l1^2;
   
   f_x = k1/(M*J) * ( -M*x(3)-m1*g*l1/2*cos(x(1)) ) - k1/M * x(3) + 1/M * ( m1*g*l1/2*(cos(x(1))*x(2)^2+sin(x(1))*x(3)) );
   g_x = k1/(M*J);
   
   u = systemController(t,x,x_ref,f_x,g_x,dof_degree);
   
   dxdt = [x(2);
           x(3);
           x(4);
           f_x+g_x*u];
end