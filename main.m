clear all
close all

%% Define the working folder
folder_name = 'data';

%% All dataset array
list=dir(folder_name);
list(1:2)=[];

%% Results (left, right, intersect distance, droplet height)
results=zeros(numel(list),4);
%% Detect bounding box
[lowx,highx,lowy,highy]=getBoundingBox(folder_name,list);
% Display bounding box
image_file = sprintf('Momentaufnahme - %02d.png',5);
rgb = imread(fullfile(folder_name,image_file));
I0 = rgb2gray(rgb);

I = I0(lowy:highy,lowx:highx);
figure
title('Cutout image')
imshow(I)
%% Droplet boundary
[droplet_boundary] = getDropletBoundaries(folder_name,list,lowy,highy,lowx,highx);
%% Surface line
close all
[p_bar] = getSurface(folder_name,list,lowy,highy,lowx,highx,droplet_boundary);
% Display surface line
figure
imshow(I,[])
hold on
plot(1:size(I,2),p_bar(2)+(1:size(I,2))*p_bar(1))
%% Surface Row
% This is for removing the depth of bar appearaing in the images. So the
% watershed algorithm works fine.
extra_row=10;
surfaceRow=round(max([p_bar(2)+1*p_bar(1) p_bar(2)+(highx-lowx)*p_bar(1)]))+extra_row;
%% Watershed transformation and fit 'very first ellipses'
% % If the result is bad, it is most likely because of the marker image.
close all
numberOfImages=size(list,1);
[ws_volume,Icutout_volume,droplet_volume,ellipse_data] = ... 
getWatershedVolume(folder_name,numberOfImages,surfaceRow,lowy,highy,lowx,highx,droplet_boundary,p_bar);
%% Get the initial ellipse as median of 'very first ellipses' for further fitting
initial_ellipse=median(ellipse_data);
ellipse_data_temp=zeros(size(list,1),5);
close all

%% Fit an ellipse again using the very first ellipse as the initial point 
% and calculate results (theta left, theta left, intersect distance, droplet height)
for i=1:size(list,1)
image_number=i;
[theta_left,theta_right,dist_droplet,dist_intersect] = getEllipse(Icutout_volume,droplet_volume,image_number,surfaceRow,p_bar,ellipse_data,initial_ellipse);
results(image_number,1)=theta_left;
results(image_number,2)=theta_right;
results(image_number,3)=dist_intersect;
results(image_number,4)=dist_droplet;
end