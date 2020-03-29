function [droplet_boundary] = getDropletBoundaries(folder_name,list,lowy,highy,lowx,highx)

numberOfImages=size(list,1);
droplet_boundaries = zeros(numberOfImages,3);
for i = 1:numberOfImages
image_file = sprintf('Momentaufnahme - %02d.png',i);
rgb = imread(fullfile(folder_name,image_file));
I = rgb2gray(rgb);
I = I(lowy:highy,lowx:highx); % Working image
I = adapthisteq(I);
I = imnlmfilt(I,'SearchWindowSize',41,'ComparisonWindowSize',11,'DegreeOfSmoothing',10);%,[15 15]);%(Ioriginal, net);

[a,b]=size(I);
Iedge = edge(I);
[row,col] = find(Iedge);
max_row = max(row);
Iedgetmp = Iedge(round(a*0.4):round(a*0.6),:);
[row,col] = find(Iedgetmp);
min_col=min(col);
max_col=max(col);
droplet_boundaries(i,:) = [min_col,max_col,max_row];
end
droplet_boundary = round(mean(rmoutliers(droplet_boundaries)));
end