function [theta_left,theta_right,dist_droplet,dist_intersect] = getEllipse(Icutout_volume,droplet_volume,image_number,surfaceRow,p_bar,ellipse_data,initial_ellipse)



Icutout=Icutout_volume(:,:,image_number);
droplet=droplet_volume(:,:,image_number);


%% Fit ellipse
centerX=ellipse_data(image_number,1);
centerY=ellipse_data(image_number,2);

major=initial_ellipse(3);
minor=initial_ellipse(4);

alpha=initial_ellipse(5);
alpha_max=max(ellipse_data(:,5)); % Maximum of all alpha values
alpha_min=min(ellipse_data(:,5)); % Minimum of all alpha values


boundary=[0.8*centerX,0.8*initial_ellipse(2),major,minor,alpha_min;
          1.2*centerX,1.2*initial_ellipse(2),major,minor,alpha_max];
      
lb=min(boundary);
ub=max(boundary);

[row,col] = find(droplet);
row=row+double(surfaceRow); % Put the bar cutout back
x = transpose([col,row]);
costfunction=@(z) Residuals_ellipse(x',z);
costIntersect = @(z) optIntersections(z,p_bar(1),p_bar(2));

out = fmincon(costfunction,[centerX,initial_ellipse(2:5)]',[],[],[],[],lb,ub,costIntersect);

z=out(1:2);
a=out(3);
b=out(4);
alpha=out(5);
m=p_bar(1);


% % To show the binarized droplet image on the cutout image
% [row,col] = find(droplet);
% idx = sub2ind(size(Icutout), row+double(surfaceRow), col);
% Icutout(idx)=255;

% % The calculation of theta angles
normal_vector=[1;m];
normal_vector=normal_vector/norm(normal_vector);
normal_vector=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*normal_vector;

c=p_bar(2);

[x1,x2,n1,n2] = getAngleBetweenLineAndEllipse(z,a,b,alpha,m,c);
% Display results
figure
imshow(Icutout,[])
hold on

plot(z(1),z(2),'ro')
plotellipse([z(1),z(2)], a, b, alpha, 'r--')
plot(1:size(Icutout,2),p_bar(2)+(1:size(Icutout,2))*p_bar(1))

t=linspace(-75,75,100);
ntmp=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*n2;
plot(x2(1)+t*ntmp(1),x2(2)+t*ntmp(2),'g')
ntmp=[cos(pi/2), -sin(pi/2);sin(pi/2), cos(pi/2)]*n1;
plot(x1(1)+t*ntmp(1),x1(2)+t*ntmp(2),'r')

theta_left=acos(dot(normal_vector,n1))/pi*180;
theta_right=acos(dot(normal_vector,n2))/pi*180;

theta_left=max(theta_left,abs(180-theta_left))
theta_right=max(theta_right,abs(180-theta_right))
%% distance between the intersection points
dist_intersect = abs(x2(1) - x1(1));
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
y_bar = x * p_bar(1) + p_bar(2);
dist_droplet = abs(y_droplet - y_bar);
plot(x1(1):x2(1), (x1(1):x2(1)) * p_bar(1) + p_bar(2),'b')
line([x x], [y_bar y_droplet]);
mkdir results;
saveas(gcf,[fullfile('results', sprintf('%03d',image_number )),'.png'])
end
