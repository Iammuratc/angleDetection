function [lowx,highx,lowy,highy]=getBoundingBox(folder_name,list)

% Define number of images to be processed to calculate the bounding box
numberOfImages = size(list,1);
lowxs=zeros(numberOfImages,1);
highxs=zeros(numberOfImages,1);
lowys=zeros(numberOfImages,1);
highys=zeros(numberOfImages,1);
    for i=1:numberOfImages
        image_number = randi([1,numel(list)]);
    image_file = sprintf('Momentaufnahme - %02d.png',image_number);
    rgb = imread(fullfile(folder_name,image_file));

    im = rgb2gray(rgb); % For display
    im=double(im);
    im=im/max(im(:));
% Cut some of the image
[b,a]=size(im);
col1=floor(0.15*a);
col2=floor(0.75*a);
row1=floor(0.12*b);
row2=floor(0.5*b);

im(:,[1:col1,col2:end])=[];

im([1:row1,row2:end],:)=[];

% Gradient image
im_smooth=imgaussfilt(im,10);
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(im_smooth, hy, 'replicate');
Ix = imfilter(im_smooth, hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

%% Find right and left cuts
beamwidth=42;
barwidth = 1300;
extra_cut = 250;
y=sum(gradmag);
% plot(y)
y=smoothdata(y,'gaussian',size(gradmag,2)/30);
[pks,locs] = findpeaks(y);
[tmp,pos]=max(pks);
beam_pos=locs(pos);
bar_right=beam_pos-beamwidth-extra_cut;
bar_left = beam_pos-beamwidth-barwidth+2*extra_cut;
% highx=beam_pos;
% im_cut=im(:,1:beam_pos);
% grad_cut=gradmag(:,1:beam_pos);
grad_cut=gradmag(:,bar_left:bar_right);
% imshow(grad_cut,[])
%% Find top and bottom buts
y=sum(grad_cut');
y=smoothdata(y,'gaussian',size(grad_cut,1)/30);
% plot(y)
[pks,locs] = findpeaks(y);
[pls,perm]=sort(pks,'descend');

max_1 = pls(1);
max_2 = pls(2);
i_max_1 = locs(perm(1));
i_max_2 = locs(perm(2));
% If there are two peaks
if abs(max_1-max_2) < 50 && abs(i_max_1-i_max_2) < 120
    if i_max_1 < i_max_2
        bar_top = i_max_2;
    else
        bar_top = i_max_1;
    end
% If there is one peak only
else
    bar_top = i_max_1;
end
 bar_bottom = bar_top + 300;
 if bar_top > 40
     bar_top = bar_top -40;
 else
     bar_top = 1;
 end
%  im_cut2=im(bar_top:bar_bottom,bar_left:bar_right);
    lowxs(i)=int16(bar_left+col1);
    highxs(i)=int16(bar_right+col1);
    lowys(i)=int16(bar_top+row1);
    highys(i)=int16(bar_bottom+row1);
 %% Display
%  figure('name',sprintf('image number:%01d', image_number));
% figure, imshow(im_cut2,[])
    end
 lowy=uint16(median(lowys));
highy=uint16(median(highys));
lowx=uint16(median(lowxs));
highx=uint16(median(highxs));
end


