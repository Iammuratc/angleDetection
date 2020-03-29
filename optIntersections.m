function [c,ceq] = optIntersections(z,m,n)
zx=z(1);
zy=z(2);
a=z(3);
b=z(4);
A=z(5);

% zx zy a b A: Ellipse coordinates
% m n: Surface line y=mx+n
% syms zx zy a b A x m n
% f = ((x-zx)*cos(A)+(m*x+n-zy)*sin(A))^2/(a)^2 + ...
%     ((x-zx)*sin(A)-(m*x+n-zy)*cos(A)).^2/(b)^2
% collect(f)
a1 = ((cos(A) + m*sin(A))^2/a^2 + (sin(A) - m*cos(A))^2/b^2);
a2 = (- (2*(cos(A) + m*sin(A))*(zx*cos(A) - sin(A)*(n - zy)))/a^2 - (2*(sin(A) - m*cos(A))*(zx*sin(A) + cos(A)*(n - zy)))/b^2);
a3 = (zx*cos(A) - sin(A)*(n - zy))^2/a^2 + (zx*sin(A) + cos(A)*(n - zy))^2/b^2 - 1;
c = -(a2^2-4*a1*a3);
ceq=[];