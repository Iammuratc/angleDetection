function [ x1, x2, y1, y2 ,n1,n2 ] = lineEllipse(a,b,m,c )

aa=(1/a^2+m^2/b^2);
bb=2*m*c/b^2;
cc=c^2/b^2-1;

x1=(-bb-sqrt(bb^2-4*aa*cc))/(2*aa);
x2=(-bb+sqrt(bb^2-4*aa*cc))/(2*aa);

y1=m*x1+c;
y2=m*x2+c;

n1=[2*x1/a^2;2*y1/b^2];
n1=n1/norm(n1);

n2=[2*x2/a^2;2*y2/b^2];
n2=n2/norm(n2);
% 
% figure
% hold on
% 
% x=linspace(0,2*pi,1000);
% plot(a*cos(x),b*sin(x))
% x=-a:a;
% plot(x,m*x+c)
% plot(x1,y1,'ro')
% plot(x2,y2,'ro')
% 
% t=linspace(0,1,2);
% plot(x1+t*40*n1(1),y1+t*40*n1(2))
% plot(x2+t*40*n2(1),y2+t*40*n2(2))
end
