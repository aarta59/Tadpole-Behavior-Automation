streletest2 = imopen(I3, strel('disk',4));
imtest3 = edge(streletest2);
figure, imshow(imtest3)

addedstel = streletest2*2;
imshow(addedstel)

%%
%FUCKING WORKING DOT DETECTION HOLY FUCK


I = imread('Frame 0100.png');

I = rgb2gray(I);
background = imopen(I, strel('disk',23));
I2 = I - background;
I3 = imadjust(I2);
streletest2 = imopen(I3, strel('disk',4));
addedstel = streletest2*2;
[centers, radius] = imfindcircles(addedstel,[8 15],'ObjectPolarity','bright', 'Sensitivity',0.90);

figure
imshow(I)
asd = viscircles(centers,radius);

%%

close all

directory = uigetdir;
cd(directory);
moviename = uigetfile('*.mov');
folder = fullfile(directory);
movFullFile = fullfile(folder, moviename);
mov = VideoReader(movFullFile);

%video dimentions
numFrames = mov.NumberOfFrames;
vidH = mov.Height;
vidW = mov.Width;

ytop = 140;
ybott = 910;
xleft = 90;
xright = 1250;

dotzeros = uint8(zeros(vidH,vidW,numFrames));

for i = 90:numFrames
    orig_img = read(mov,i);
    orig_img = rgb2gray(orig_img);
    croped_orig = orig_img(ytop:ybott,xleft:xright,:);
    
    backg = imopen(croped_orig, strel('disk',23));
    minus_bck = croped_orig - backg;
    adj_mius = imadjust(minus_bck);
    str_mius = imopen(adj_mius, strel('disk',4));
    str_mius = str_mius*2;
    dotzeros(ytop:ybott,xleft:xright,i) = str_mius;
end

Q_loc_estimateX(isnan(Q_loc_estimateX)) = [];
Q_loc_estimateY(isnan(Q_loc_estimateY)) = [];

allcenter = cell(1,numFrames);
allradius = cell(1,numFrames);
for i = 90:numFrames
    [allcenter{i}, allradius{i}] = imfindcircles(dotzeros(:,:,i),[8 15],...
        'ObjectPolarity','bright', 'Sensitivity',0.90);    
end
 
for i = 90:length(Q_loc_estimateX)
%     mov_img = read(mov,i+14);
%     mov_img = rgb2gray(mov_img);
%     imshow(mov_img)
    hold on
    plot(Q_loc_estimateY(i),Q_loc_estimateX(i),'og')
    d = allradius{i+14}*2;
    px = allcenter{i+14}(:,1) - allradius{i+14};
    py = allcenter{i+14}(:,2) - allradius{i+14};
    for j = 1:length(d)
        h = rectangle('Position',[px(j) py(j) d(j) d(j)],'Curvature',[1,1],'FaceColor',[0,0,0]);
        
    end
    set(gca,'Ydir','reverse')
    axis([0 1344 0 1024])
    pause(0.1)
    clf
end
 



save('dot_centers_radii.mat','allcenter','allradius')    
    
    
load('dot_centers_radii.mat')
    
    
    