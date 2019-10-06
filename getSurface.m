function [p_left,p_right,min_col,mid,max_col,max_row] = getSurface(folder_name,v_folder,image_number,row0,row1,col0,col1)

image_file = sprintf('Momentaufnahme - %02d.png',image_number);
rgb = imread(fullfile('images',folder_name,v_folder,image_file));
Iextended = rgb2gray(rgb); % For display
Icutout = Iextended(row0:row1,col0:col1); % Working image

%net = denoisingNetwork('DnCNN');

Icutout = adapthisteq(Icutout);
Icutout = imnlmfilt(Icutout,'SearchWindowSize',41,'ComparisonWindowSize',11,'DegreeOfSmoothing',10);%,[15 15]);%(Ioriginal, net);
I = Icutout;
[a,b]=size(I);

%% Droplet boundaries
Iedge = edge(I);
% figure,imshow(Iedge);

[row,col] = find(Iedge);
max_row = max(row);
Iedgetmp = Iedge(round(a*0.4):round(a*0.6),:);
[row,col] = find(Iedgetmp);
min_col=min(col);
max_col=max(col);
mid=round((min_col+max_col)/2);
%% Surface
Itemp=imbinarize(I);
% figure,imshow(Itemp)
Ileft=I(:,1:min_col-10);
Iright=I(:,max_col+10:end);
Ileft=imbinarize(Ileft);
Iright=imbinarize(Iright);
Ileft=~Ileft;
Iright=~Iright;
Iright=imdilate(Iright,strel('disk',1))-Iright;
Ileft=imdilate(Ileft,strel('disk',1))-Ileft;
% figure,imshow(Ileft)
% figure,imshow(Iright)
[row_l,col_l] = find(Ileft);
[row_r,col_r] = find(Iright);
% 
p_left = polyfit(col_l,row_l,1);
p_right = polyfit(col_r,row_r,1);
end

