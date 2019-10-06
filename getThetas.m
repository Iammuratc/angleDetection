function [theta_left,theta_right,dist_droplet,dist_intersect]=getThetas(droplet,Icutout,p_left,p_right,surfaceRow)
% Iellipse=zeros(size(Icutout));
% [row,col] = find(Ileft);
% idx = sub2ind(size(Icutout), row+surfaceRow, col);
% % Icutout(idx)=255;
% Iellipse(idx)=255;
% [row,col] = find(Iright);
% idx = sub2ind(size(Icutout), row+surfaceRow, col);
% % Icutout(idx)=255;
% Iellipse(idx)=255;
% figure,imshow(Iellipse)

%% Boundary conditions
Iellipse=droplet;
[row,col] = find(Iellipse);
row=row+double(surfaceRow);
x0=col;
y0=row;
z0=mean([x0,y0]);
a0=max(x0)-min(x0);
b0=max(y0)-min(y0);
a0=a0/1.6;
b0=b0/2.2;
lb=[min(x0), min(y0),0.75*a0,0.75*b0,-pi];
hb=[max(x0),max(y0),1.5*a0,1.5*b0,pi];
% For surface line draw
mi=min(x0);
ma=max(x0);
mid=round((mi+ma)/2);

%% Fit one ellipse only
x = transpose([col,row]);
costfunction=@(z) Residuals_ellipse(x',z);
costIntersect = @(z) optIntersections(z,p_left(1),p_left(2));
% options = optimoptions('fmincon','Display','iter');
out = fmincon(costfunction,[z0(1),z0(2),a0,b0,0]',[],[],[],[],lb,hb,costIntersect);

z=out(1:2);
a=out(3);
b=out(4);
alpha=out(5);
% figure
% imshow(Icutout,[])
% hold on
% axis equal
% plotellipse(z, a, b, alpha, 'r--')
%% Theta angles
% 
m_left=p_left(1);
m_right=p_right(1);
% % 
normal_vector_left=[1;m_left];
normal_vector_left=normal_vector_left/norm(normal_vector_left);

normal_vector_left=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*normal_vector_left;


normal_vector_right=[1;m_right];
normal_vector_right=normal_vector_right/norm(normal_vector_right);

normal_vector_right=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*normal_vector_right;

c_left=p_left(2);
c_right=p_right(2);

[x1_left,x2_left,n1_left,n2_left] = getAngleBetweenLineAndEllipse(z,a,b,alpha,m_left,c_left);
[x1_right,x2_right,n1_right,n2_right] = getAngleBetweenLineAndEllipse(z,a,b,alpha,m_right,c_right);
x1=x1_left;
n1=n1_left;

x2=x2_right;
n2=n2_right;

%% Display results
figure
imshow(Icutout,[])
hold on

plotellipse([z(1),z(2)], a, b, alpha, 'r--')
% plotellipse(z_right, a_right, b_right, alpha_right, 'g--')
plot(mi:mid,p_left(2)+[mi:mid]*p_left(1))
plot(mid:ma,p_right(2)+[mid:ma]*p_right(1))
% % 
plot(x1_left(1),x1_left(2),'ro')
plot(x2_left(1),x2_left(2),'ro')
plot(x1_right(1),x1_right(2),'go')
plot(x2_right(1),x2_right(2),'go')
t=linspace(-75,75,100);
ntmp=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*n2;
plot(x2(1)+t*ntmp(1),x2(2)+t*ntmp(2),'g')
ntmp=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*n1;
plot(x1(1)+t*ntmp(1),x1(2)+t*ntmp(2),'r')
% 
% 
% 
% % 
theta_left=acos(dot(normal_vector_left,n1))/pi*180;
theta_right=acos(dot(normal_vector_right,n2))/pi*180;
% % 
theta_left=max(theta_left,abs(180-theta_left))
theta_right=max(theta_right,abs(180-theta_right))
%% distance between the intersection points
dist_intersect = abs(x2_right(1) - x1_left(1));
%% distance between the droplet top point and the bar surface
% save(sprintf('ellipse_data_%02d',image_number),'z', 'a', 'b', 'alpha')
npts = 100;
t = linspace(0, 2*pi, npts);

% Rotation matrix
Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
% Ellipse points
X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);
y = X(2,:);
[y_droplet,pos]=max(y);
x = X(1,pos);
y_bar = x * p_left(1) + p_left(2);
dist_droplet = abs(y_droplet - y_bar);
plot(x1_left(1):x2_right(1), (x1_left(1):x2_right(1)) * p_left(1) + p_left(2),'b')
line([x x], [y_bar y_droplet]);
end


% 
% % z(2) = z(2);
% idx = sub2ind(size(Icutout), row, col);
% Icutout(idx)=255;
% figure
% imshow(Icutout,[])
% hold on
% axis equal
% plotellipse(z, a, b, alpha, 'r--')
% plot(mi:mid,p_left(2)+(mi:mid)*p_left(1))
% plot(mid:ma,p_right(2)+(mid:ma)*p_right(1))