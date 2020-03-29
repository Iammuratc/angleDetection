function [ws_volume,Icutout_volume,droplet_volume,ellipse_data] = getWatershedVolume(folder_name,numberOfImages,surfaceRow,lowy,highy,lowx,highx,droplet_boundary,p_bar)

row0=lowy;
row1=highy;
col0=lowx;
col1=highx;
max_col = droplet_boundary(2);
min_col = droplet_boundary(1);
max_row = droplet_boundary(3);

ws_volume=zeros(row1-surfaceRow-row0+1, max_col+51,numberOfImages);
% Iwork_volume=zeros(row1-surfaceRow-row0+1, max_col+51,numberOfImages);
Icutout_volume=zeros(row1-row0+1, col1-col0+1,numberOfImages);
droplet_volume=zeros(71, col1-col0+1,numberOfImages);
ellipse_data=zeros(numberOfImages,5);
for image_number =1:numberOfImages
    
image_file = sprintf('Momentaufnahme - %02d.png',image_number);
rgb = imread(fullfile(folder_name,image_file));
Ioriginal = rgb2gray(rgb); % For display
Icutout = Ioriginal(row0:row1,col0:col1);
%% Working image (left,right)
min_col = droplet_boundary(1);
max_col = droplet_boundary(2);
% extra=10;

I = Ioriginal(row0+surfaceRow:row0+surfaceRow+70,col0:col1);

I = adapthisteq(I);
I = imnlmfilt(I,'SearchWindowSize',41,'ComparisonWindowSize',11,'DegreeOfSmoothing',10);%,[15 15]);%(Ioriginal, net);

%% Watershed images
g1=[5,5,5; -3,0,-3; -3,-3,-3];
g2=[5,5,-3; 5,0,-3; -3,-3,-3];
g3=[5,-3,-3; 5,0,-3; 5,-3,-3];
g4=[-3,-3,-3; 5,0,-3; 5,5,-3];
g5=[-3,-3,-3; -3,0,-3; 5,5,5];
g6=[-3,-3,-3; -3,0,5;-3,5,5];
g7=[-3,-3,5; -3,0,5;-3,-3,5];
g8=[-3,5,5; -3,0,5;-3,-3,-3];
x1=imfilter(I,g1,'replicate','same');
x2=imfilter(I,g2,'replicate','same');
x3=imfilter(I,g3,'replicate','same');
x4=imfilter(I,g4,'replicate','same');
x5=imfilter(I,g5,'replicate','same');
x6=imfilter(I,g6,'replicate','same');
x7=imfilter(I,g7,'replicate','same');
x8=imfilter(I,g8,'replicate','same');


ws_image_left = abs(x1)+abs(x2)+abs(x3)+abs(x4)+abs(x5);

ws_image_right = abs(x5)+abs(x6)+abs(x7)+abs(x8);
%% Markers (left,right)
Iedge = edge(I,'Prewitt');

Iedge = bwareaopen(Iedge,20);
% [row_edge,col_edge]=find(Iedge);

% sorted_col = sort(col_edge);
% max_edge = sorted_col(end-30);

% Markers (right)
[row,col]=size(I);
marker_right = zeros(size(I));

% Background
marker_right(:,round(col*0.95):col) = 1;
% Droplet
marker_right(round(row*0.2):end,1:max_col-60) = 1;

% Markers (left)
marker_left = zeros(size(I));
% Background
marker_left(:,1:round(col*0.05)) = 1;
% Droplet
marker_left(:,min_col+70:end) = 1;
%% Watershed
Kleft = imimposemin(ws_image_left,marker_left);

Kright = imimposemin(ws_image_right,marker_right);

Lleft = watershed(Kleft);
Lright = watershed(Kright);


%% Display watershed results
% figure,imshow(Kleft,[])
% figure,imshow(Kright,[])
%  Lrgb = label2rgb(Lright, 'jet', 'w', 'shuffle');
% figure
% imshow(I)
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.2;
% title('Lright')
%  Lrgb = label2rgb(Lleft, 'jet', 'w', 'shuffle');
%  figure
% imshow(I)
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.2;
% title('Lleft')
%% Droplet boundaries
Ileft = Lleft==Lleft(1,1);
Ileft=Ileft-imerode(Ileft,strel('disk',1));

Iright = Lright==Lright(1,1);
Iright=Iright-imerode(Iright,strel('disk',1));


% figure, imshow(Ileft)
% figure, imshow(Iright)

droplet=zeros(size(I));
Ileft=Ileft>0;
droplet(Ileft)=1;
Iright=Iright>0;
droplet(Iright)=1;

% figure,imshow(droplet)
droplet_volume(:,:,image_number)=droplet;
%% Fit ellipse to outer boundaries
[row,col] = find(droplet);
row=row+double(surfaceRow);

r = round((max_col-min_col)/2);
% x0=mean(col);
% y0=mean(row);
centerX = mean(col);%min_col+r+5;
centerY = max(row);%max_row-r-double(surfaceRow);

lb=[0.75*centerX, 0.75*centerY,0.75*r,0.75*r,-pi];
ub=[1.5*centerX,1.5*centerY,1.5*r,1.5*r,pi];



x = transpose([col,row]);
costfunction=@(z) Residuals_ellipse(x',z);
costIntersect = @(z) optIntersections(z,p_bar(1),p_bar(2));



out = fmincon(costfunction,[centerX,centerY,r-5,r-5,0]',[0 0 1 -1 0],0,[],[],lb,ub,costIntersect);

z=out(1:2);
a=out(3);
b=out(4);
alpha=out(5);
m=p_bar(1);


% [row,col] = find(droplet);
% idx = sub2ind(size(Icutout), row+double(surfaceRow), col);
% Icutout(idx)=255;

ellipse_data(image_number,:)=[z(1) z(2) a b alpha];
%% Display results
figure
imshow(Icutout,[])
hold on

plot(centerX,centerY,'ro')
% plotellipse([centerX,centerY+double(surfaceRow)],r,r,0,'g--') %% initial circle
plotellipse([z(1),z(2)], a, b, alpha, 'r--')
plot(1:size(Icutout,2),p_bar(2)+(1:size(Icutout,2))*p_bar(1))

Icutout_volume(:,:,image_number)=Icutout;
end

% ws_left=cat(3,ws_volume(:,:,4),ws_volume(:,:,3),ws_volume(:,:,2));
% ws_right=cat(3,ws_volume(:,:,49),ws_volume(:,:,48),ws_volume(:,:,47));
% 
% ws_volume=cat(3,ws_left,ws_volume,ws_right);
% 
% sigma=2;
% n=3;
% 
% 
% filter=-n:n;
% filter=1/sqrt(2*pi*sigma^2)*exp(-filter.^2/(2*sigma^2));
% filter=filter/sum(filter(:));
% 
% ws_volume=permute(ws_volume,[1 3 2]);
% 
% ws_volume_fl=convn(ws_volume,filter,'valid');
% ws_volume=permute(ws_volume_fl,[1 3 2]);


%%
% filter=zeros(1,1,2*n+1);
% filter(1,1,:)=-n:n;
% sigma=2;
% n=3;
% filter=-n:n;
% filter=1/sqrt(2*pi*sigma^2)*exp(-filter.^2/(2*sigma^2));
% filter=filter/sum(filter(:));
% % 
% ws_volume_new=zeros(size(ws_volume));
% for k=1:50
%     for l=-n:n
%         ind=k+l;
%         if(ind<1)
%             ind=1;
%         elseif(ind>50)
%             ind=50-(ind-50);
%         end
%         ws_volume_new(:,:,k)=ws_volume_new(:,:,k)+filter(l+n+1)*ws_volume(:,:,ind);
%         
%     end
% end
% ws_volume=ws_volume_new;

