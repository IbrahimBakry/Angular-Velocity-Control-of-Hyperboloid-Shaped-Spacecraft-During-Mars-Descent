function []=DE_ALL_Forms_RK4(tspan,y0)

h = (tspan(2)-tspan(1))/100; y1(1,:) = y0;
i = 0;
for x = tspan(1):h:tspan(2)
    i = i + 1;
    % RK4 - First Form
    k11 = Fun(x,y1(i,:));
    k21 = Fun(x+0.5*h,y1(i,:)+0.5*h*k11);
    k31 = Fun(x+0.5*h,y1(i,:)+0.5*h*k21);
    k41 = Fun(x+h,y1(i,:)+h*k31);
    y1(i+1,:) = y1(i,:) + (1/6)*h*(k11 + 2*k21 + 2*k31 + k41);
    xx(i) = x;
end
end