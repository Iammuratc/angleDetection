function [droplet3,Icutout3,r_init_options,K] = getDroplet3D(folder_name,v_folder,numberOfImages,surfaceRow,lowy,highy,lowx,highx,droplet_boundary,p_bar)



% numberOfImages=1;

row0=lowy;
row1=highy;
col0=lowx;
col1=highx;
max_col = droplet_boundary(2);
min_col = droplet_boundary(1);
max_row = droplet_boundary(3);

ws_image=zeros(row1-row0-surfaceRow+1,max_col+50+1,numberOfImages);
marker=zeros(row1-row0-surfaceRow+1,max_col+50+1,numberOfImages);
Icutout3=zeros(row1-row0+1,col1-col0+1,numberOfImages);


for image_number=1:numberOfImages
% my_im = randi([1,numel(list)]);
image_file = sprintf('Momentaufnahme - %02d.png',image_number);
rgb = imread(fullfile('images',folder_name,v_folder,image_file));
Ioriginal = rgb2gray(rgb); % For display
Icutout = Ioriginal(row0:row1,col0:col1);
%% Working image
Iwork = Ioriginal(row0+surfaceRow:row1,col0:col0+max_col+50);
Iwork = adapthisteq(Iwork);
Iwork = imnlmfilt(Iwork,'SearchWindowSize',41,'ComparisonWindowSize',11,'DegreeOfSmoothing',10);%,[15 15]);%(Ioriginal, net);
%% Outer droplet boundaries
Iedge = edge(Iwork);
% figure,imshow(Ibinary)

[Idt,IDX] = bwdist(~Iedge,'euclidean');
J = Idt;
% figure,imshow(J)
[row,col]=find(J);
droplet=zeros(size(J));
for i=1:max(row)
    my_array = J(i,:);
    [row,col]=find(my_array);
    % TODO
    col_r = max(col);
    col_l = min(col);
    droplet(i,col_r) =1;
    droplet(i,col_l) =1;
end
% figure,imshow(droplet)

%% Boundary conditions
[ad,bd] = size(droplet);
droplettmp=zeros(size(droplet));
droplettmp(round(ad*0.2):round(ad*0.6),:)=droplet(round(ad*0.2):round(ad*0.6),:);
[row,col] = find(droplet);
[row_x,col_x] = find(droplettmp);
x0=col_x;
y0=row;
z0=mean([col,row]);
% a0=max(x0)-min(x0);
% b0=max(y0)-min(y0);
% a0=a0/2;
% b0=b0/2;
% lb=[min(x0), min(y0),0.75*a0,0.75*b0,-pi];
% hb=[max(x0),max(y0),1.5*a0,1.5*b0,pi];

%% Fit ellipse to outer boundaries
[imageSizeY,imageSizeX]=size(droplet);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
centerX = min_col + round((max(x0)-min(x0))/2);
centerY = round(max(y0)/2);
r=max(round((max(x0)-min(x0))/2),(max_row/2));
radius = r+50;
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
droplet(~circlePixels) =0;

lb=[min(x0), min(y0),0.75*r,0.75*r,-pi];
hb=[max(x0),max(y0),1.5*r,1.5*r,pi];

% figure,imshow(droplet)
[row,col] = find(droplet);
x = transpose([col,row]);
costfunction=@(z) Residuals_ellipse(x',z);
costIntersect = @(z) optIntersections(z,p_bar(1),p_bar(2));

out = fmincon(costfunction,[centerX,centerY,r-5,r-5,0]',[],[],[],[],lb,hb,costIntersect);
    
z=out(1:2);
a=out(3);
b=out(4);
alpha=out(5);
r_init_options=[a b r];
%% Display
% figure
% imshow(droplet,[])
% hold on
% axis equal
% plotellipse(z, a, b, alpha, 'r--')
%% Markers
[imageSizeY,imageSizeX]=size(droplet);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
centerX = z(1);
centerY = z(2);
% Droplet
r=mean([a,b]);
radius = r-10;
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;

Imarker=marker(:,:,image_number);

[row,col] = find(circlePixels);
idx = sub2ind(size(Imarker), row, col);
Imarker(idx)=255;
marker(:,:,image_number)=Imarker;

% Background
radius = r+25;
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;

[row,col] = find(~circlePixels);
idx = sub2ind(size(Imarker), row, col);
Imarker(idx)=255;
marker(:,:,image_number)=Imarker;

%% Watershed image
g1=[5,5,5; -3,0,-3; -3,-3,-3];
g2=[5,5,-3; 5,0,-3; -3,-3,-3];
g3=[5,-3,-3; 5,0,-3; 5,-3,-3];
g4=[-3,-3,-3; 5,0,-3; 5,5,-3];
g5=[-3,-3,-3; -3,0,-3; 5,5,5];
g6=[-3,-3,-3; -3,0,5;-3,5,5];
g7=[-3,-3,5; -3,0,5;-3,-3,5];
g8=[-3,5,5; -3,0,5;-3,-3,-3];

x1=imfilter(Iwork,g1,'replicate','same');
x2=imfilter(Iwork,g2,'replicate','same');
x3=imfilter(Iwork,g3,'replicate','same');
x4=imfilter(Iwork,g4,'replicate','same');
x5=imfilter(Iwork,g5,'replicate','same');
x6=imfilter(Iwork,g6,'replicate','same');
x7=imfilter(Iwork,g7,'replicate','same');
x8=imfilter(Iwork,g8,'replicate','same');

I = abs(x8)+abs(x7)+abs(x2)+abs(x1);
ws_image(:,:,image_number)=I;

Icutout3(:,:,image_number) = Icutout;
% figure,imshow(imimposemin(I,marker(:,:,image_number)),[])
end
%% Watershed
K= imimposemin(ws_image,marker);
L = watershed(K);

% figure,imshow(K,[])
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure
% imshow3(Icutout3)
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.2;
% title('Lrgb superimposed transparently on original image')
%% Fit ellipse
droplet3 = L==L(round(centerY),round(centerX));
droplet3=droplet3-imerode(droplet3,strel('disk',1));
% Display
% figure,imshow(droplet3(:,:,5),[])
%% Display on the original image

% ind = find(droplet3);
% [idx1,idx2,idx3] = sub2ind(size(Icutout3), ind);%
% Icutout3(idx1,idx2,idx3)=255;
% figure,imshow3(Icutout3)
%%
