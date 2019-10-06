function [droplet,Icutout] = getDroplet(folder_name,v_folder,image_number,surfaceRow,lowy,highy,lowx,highx,droplet_boundary,p_left)

row0=lowy;
row1=highy;
col0=lowx;
col1=highx;
max_col = droplet_boundary(2);

image_file = sprintf('Momentaufnahme - %02d.png',image_number);
rgb = imread(fullfile('images',folder_name,v_folder,image_file));
Ioriginal = rgb2gray(rgb); % For display
Icutout = Ioriginal(row0:row1,col0:col1);
%% Working image
Iwork = Ioriginal(row0+surfaceRow:row1,col0:col0+max_col+20);
figure,imshow(Icutout)
Iwork = adapthisteq(Iwork);
Iwork = imnlmfilt(Iwork,'SearchWindowSize',41,'ComparisonWindowSize',11,'DegreeOfSmoothing',10);%,[15 15]);%(Ioriginal, net);
%% Outer droplet boundaries
Ibinary = edge(Iwork);
% figure,imshow(Ibinary)

[Idt,IDX] = bwdist(~Ibinary,'euclidean');
J = Idt;
[row,col]=find(J);
droplet=zeros(size(J));
% row = unique(row);
for i=1:max(row)
    my_array = J(i,:);
    [row,col]=find(my_array);
    col_r = max(col);
    col_l = min(col);
    droplet(i,col_r) =1;
    droplet(i,col_l) =1;
end
% figure,imshow(droplet)
%% Boundary conditions
[row,col] = find(droplet);
x0=col;
y0=row;
z0=mean([x0,y0]);
a0=max(x0)-min(x0);
b0=max(y0)-min(y0);
a0=a0/3;
b0=b0/3;
lb=[min(x0), min(y0),0.75*a0,0.75*b0,-pi];
hb=[max(x0),max(y0),1.5*a0,1.5*b0,pi];
% lb=[min(x0), min(y0),0.5*a0,0.5*b0,-pi];
% hb=[max(x0),max(y0),2*a0,2*b0,pi];
% For surface line draw
mi=min(x0);
ma=max(x0);
mid=round((mi+ma)/2);

%% Fit ellipse to outer boundaries
x = transpose([col,row]);
costfunction=@(z) Residuals_ellipse(x',z);
costIntersect = @(z) optIntersections(z,p_left(1),p_left(2));
% options = optimoptions('fmincon','Display','iter');
try
    out = fmincon(costfunction,[z0(1),z0(2),a0,b0,0]',[],[],[],[],lb,hb,costIntersect);
catch

    out = fmincon(costfunction,[z0(1),z0(2),a0,b0,0]',[],[],[],[],lb,hb,costIntersect);
end
    
z=out(1:2);
a=out(3);
b=out(4);
alpha=out(5);
% Display
% figure
% imshow(Iwork,[])
% hold on
% axis equal
% plotellipse(z, a, b, alpha, 'r--')
%% Markers (Circle)
[imageSizeY,imageSizeX]=size(droplet);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
centerX = z(1);
centerY = z(2);
% Droplet
radius = a-7;
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
marker =zeros(size(droplet));
marker(circlePixels) =1;
% Background
radius = a+5;
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
marker(~circlePixels) =1;

% Display
% image(circlePixels) ;
% colormap([0 0 0; 1 1 1]);
% title('Binary image of a circle');

%% Watershed image
Itemp=imbinarize(Iwork);
[Idt,IDX] = bwdist(~Itemp,'euclidean');
ws_image = Idt;


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

I = abs(x1)+abs(x2)+abs(x7)+abs(x8);%+abs(x4)+abs(x5)+abs(x6);%abs(x1)+abs(x2)+abs(x7)+abs(x8)
ws_image=I;
% figure,imshow(ws_image,[])
%% Watershed
K= imimposemin(ws_image,marker);
% figure,imshow(K)
L = watershed(K);
% Display
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure
% imshow(Iwork)
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.1;
% title('L')
%% Fit ellipse
droplet = L==L(round(centerY),round(centerX));
droplet=droplet-imerode(droplet,strel('disk',1));
% Display
% figure,imshow(droplet)