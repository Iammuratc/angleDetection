function [x1,x2,n1,n2] = getAngleBetweenLineAndEllipse(z,a,b,theta,m,c)


p0=[0;c];
n=[1;m];
n=n/norm(n);
n=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*n;

p0=p0-z;
R=[cos(theta), -sin(theta);sin(theta), cos(theta)]';
n=R*n;
p0=R*p0;

m_2=-n(1)/n(2);
c_2=n(1)*p0(1)/n(2)+p0(2);



[ x1, x2, y1, y2,n1,n2 ] = lineEllipse(a,b,m_2,c_2 );

% 
% figure
% hold on
% 
% x=linspace(0,2*pi,1000);
% plot(a*cos(x),b*sin(x))
% x=-a:a;
% plot(x,m_2*x+c_2)
% plot(x1,y1,'ro')
% plot(x2,y2,'ro')
% t=linspace(0,1,2);
% plot(x1+t*40*n1(1),y1+t*40*n1(2))
% plot(x2+t*40*n2(1),y2+t*40*n2(2))
% daspect([1 1 1])


x1=[x1;y1];
x2=[x2;y2];

x1=R'*x1+z;
x2=R'*x2+z;
n1=R'*n1;
n2=R'*n2;

% figure
% x=linspace(0,2*pi,100);
% plot(a*cos(x),b*sin(x));
% hold on
% plot(-a:a, c_2+m_2*[-a:a])

% figure
% plotellipse(z,a,b,theta);
%  hold on
% 
% % 
% %x1
% %x2
% t=linspace(0,1,2);
% 
%  %plot(-a:a, c+m*[-a:a])
% plot(x1(1)+t*40*n1(1),x1(2)+t*40*n1(2))
% plot(x2(1)+t*40*n2(1),x2(2)+t*40*n2(2))
% daspect([1 1 1])

end

