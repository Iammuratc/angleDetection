clear all
close all

%% Define the working folder

folder_name = 'Sphärischer Graphit (Sorte SG)';
v_folder = 'KWSG_V011_MA_2005';

list=dir(fullfile('images',folder_name,v_folder));
list(1:3)=[];

%% Detect bounding box
numberOfImages = 5;
[lowx,highx,lowy,highy]=getBoundingBox(folder_name,v_folder,list,numberOfImages);
% Display
image_file = sprintf('Momentaufnahme - %02d.png',2);
rgb = imread(fullfile('images',folder_name,v_folder,image_file));
I = rgb2gray(rgb); 
I = I(lowy:highy,lowx:highx);
figure, imshow(I)
%% Detect the average line
numberOfImages = 10;
results = zeros(numel(list),4);

p_lefts = zeros(numberOfImages,2);
p_rights = zeros(numberOfImages,2);
mi_list = zeros(numberOfImages,1);
mid_list = zeros(numberOfImages,1);
ma_list = zeros(numberOfImages,1);
droplet_boundaries = zeros(numberOfImages,3);
for i = 1:numberOfImages
    image_number = randi([1,numel(list)]);
    [p_left,p_right,min_col,mid,max_col,max_row] = getSurface(folder_name,v_folder,image_number,lowy,highy,lowx,highx); 
    droplet_boundaries(i,:) = [min_col,max_col,max_row];
    p_lefts(i,:) = p_left;
    p_rights(i,:) = p_right;
    mi_list(i) = min_col;
    mid_list(i) = mid;
    ma_list(i) = max_col;
end
p_lefts = rmoutliers(p_lefts);
p_rights = rmoutliers(p_rights);
p_left = mean(p_lefts);
p_right = mean(p_rights);
% p_average_left = (p_average_left + p_average_right)/2; % Get average for left and right lines (delete if you dont want)
% p_average_right = p_average_left; % Get average for left and right lines (delete if you dont want)

min_col = mean(mi_list);
mid = mean(mid_list);
max_col = mean(ma_list);
droplet_boundary = round(mean(rmoutliers(droplet_boundaries)));
%% Plot surface line
figure
imshow(I,[])
hold on
plot(min_col:mid,p_left(2)+[min_col:mid]*p_left(1))
plot(mid:max_col,p_right(2)+[mid:max_col]*p_right(1))
surfaceRow = round(max([p_left(2)+1*p_left(1) ...
    p_right(2)+(highx-lowx)*p_right(1)]))+5;
%% Droplet detection
close all

for image_number =1:numel(list)
    disp(image_number)
    [droplet,Icutout] = getDroplet(folder_name,v_folder,image_number,surfaceRow,lowy,highy,lowx,highx,droplet_boundary,p_left);
    [theta_left,theta_right,dist_droplet,dist_intersect]=getThetas(droplet,Icutout,p_left,p_right,surfaceRow);
    results(image_number,1)=theta_left;
    results(image_number,2)=theta_right;
    results(image_number,3)=dist_intersect;
    results(image_number,4)=dist_droplet;
end