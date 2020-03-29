function [p_surface] = getSurface(folder_name,list,lowy,highy,lowx,highx,droplet_boundary)

max_col=droplet_boundary(2);
min_col = droplet_boundary(1);

% Define number of images to be processed to calculate the surface line
numberOfImages=4;
Isum=uint16(zeros(71,max_col+151));

for i =1:numberOfImages
image_number = randi([1,numel(list)]);
image_file = sprintf('Momentaufnahme - %02d.png',image_number);
rgb = imread(fullfile(folder_name,image_file));
I = rgb2gray(rgb);
I = I(lowy:lowy+70,lowx-150:lowx+max_col);

%% Contrast enhancement
I = adapthisteq(I);
I = imnlmfilt(I,'SearchWindowSize',41,'ComparisonWindowSize',11,'DegreeOfSmoothing',10);%,[15 15]);%(Ioriginal, net);


%% Add up all images
Isum(:,:)=Isum+uint16(I);
figure,imshow(edge(I))
end
Isum = imadjust(Isum,stretchlim(Isum),[]);


%% Get the edge image
Iedge = edge(Isum,'Prewitt');
%% Hough transform
[H,T,R] = hough(Iedge,'Theta',-90:0.01:-89);

P  = houghpeaks(H,1);

x=linspace(1,size(I,2),1000);
my_rho=R(P(:,1));
my_theta=deg2rad(T(P(:,2)));
y=(my_rho-x*cos(my_theta))/sin(my_theta);

% figure
% imshow(Iedge)
% hold on
% plot(x,y,'r')

p_surface = [cos(my_theta)/sin(my_theta), my_rho/sin(my_theta)];
end
