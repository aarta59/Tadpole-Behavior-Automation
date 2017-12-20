%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Tadpole Tracker 9000                           %
%                 Christopher Marshall, November 2017                   %
%            For Cline Lab, Dorris Center for Neuroscience              %
%           Scripps Research Institute, La Jolla, California            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
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

%Initializations 
fullVideo = zeros(vidH,vidW,numFrames);
img = zeros(vidH,vidW,numFrames);
noDot_img = zeros(vidH,vidW,numFrames);

%Reading and saving all frames of video to 3d matrix
for i = 1:numFrames
    vid_img = read(mov,i);
    fullVideo(:,:,i) = vid_img(:,:,1);
end

%Taking Median of video frames
for i = 1:numFrames
    img_tmp = read(mov,i);
    img_tmp = rgb2gray(img_tmp);
    img_tmp = imopen(img_tmp, strel('disk',25));
    img(:,:,i) = img_tmp;
end

bck_img = (mean(img,3));
bck_img = uint8(bck_img);

%Removing dots from each frame and saving to new matrix
for i = 1:numFrames
    orig_img = read(mov,i);
    orig_img = rgb2gray(orig_img);
    dot_str = orig_img - bck_img;
    noDot_tmp = orig_img - dot_str;
    noDot_img(:,:,i) = noDot_tmp;
end

% Initialize log gaussian filter
%for example video 95 hsizeh=60 and sigmah=8

hsizeh = 60;  
sigmah = 8;   
h = fspecial('log', hsizeh, sigmah);

clear i
%% Iteratively finding tadpoles from blobs

X = cell(1,numFrames); %detection X coordinate indice
Y = cell(1,numFrames);  %detection Y coordinate indice

for i = 15:numFrames
    %img_real = (read(mov,i));
    bck_img = double(bck_img);
    img = noDot_img(:,:,i); 
    sub_img = (img - bck_img);

    %Blob Filtering
    blob_img = filter2(h,sub_img,'same');
    
    %Thresholding level for blob
    idx = find(blob_img < 0.03); 
    blob_img(idx) = nan;
    
    %http://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m-extrema2-m
    %Finds peak indices for blobs
    [zmax,imax,zmin,imin] = extrema2(blob_img);
    [X{i},Y{i}] = ind2sub(size(blob_img),imax);

    %Plot of raw detections with threshold overlay
%     imagesc(blob_img)
%     hold on
%     for j = 1:length(X{i})
%        plot(Y{i}(j),X{i}(j),'or')
%     end
%     axis off
%     pause
     
end

save('raw_tad_detections.mat','X','Y')

%% Clear Memory Before Assignment Begins



%% Loading Detections Matrix

%load('raw_tad_detections.mat')

%% Kalman Filter Variable Definitions 

dt = 1; %sampling rate
s_frame = 15; %detection starting frame


u = 0; %starting acceleration magnitude 
Tad_noise_mag = 1; %variability in tadpole speed [=]m/s^2
tmn_x = 0.1; %noise in horizontal direction, x-axis
tmn_y = 0.1; %noise in vertical direction, y-axis

%Process noise into covariance matrix (Ex)
Ez = [tmn_x 0; 0 tmn_y];
Ex = [dt^4/4 0 dt^3/2 0; 0 dt^4/4 0 dt^3/2; dt^3/2 0 dt^2 0; 0 dt^3/2 0 dt^2].*Tad_noise_mag^2;

%Estimate of initial position variance (covariance matrix)
P = Ex;

%2-D updates for coefficent matrix
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrice
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0];


%Initializing results
Q_loc_measure = [];

%Initialize estimations in 2-D
Q = [X{s_frame} Y{s_frame} zeros(length(X{s_frame}),1) zeros(length(X{s_frame}),1)]';
Q_est = nan(4,2000);
Q_est(:,1:size(Q,2)) = Q; %initial location estimate
Q_loc_estimateY = nan(2000);
Q_loc_estimateX = nan(2000);
P_estimate = P;
numDet = size(X{s_frame},1); %number of detections
numF = find(isnan(Q_est(1,:))==1, 1)-1; %number of estimates


for t = s_frame:(numFrames)-1
    
    %load detections matrix
    Q_loc_measure = [X{t} Y{t}];
    
    %Kalman Filter
    numDet = size(X{t},1);
    for F = 1:numF
        Q_est(:,F) = A * Q_est(:,F) + B * u;
    end
    
    %Predict covariance
    P = A * P * A' + Ex;
    
    %Kalman gain
    K = P * C' * inv(C * P * C' + Ez);
    
    %Assign detections (Hungarian/Munkres Algorithm)
    est_dist = pdist([Q_est(1:2,1:numF)'; Q_loc_measure]);
    est_dist = squareform(est_dist);
    est_dist = est_dist(1:numF, numF+1:end);
    
    [asign, cost] = munkres(est_dist);
    asign = asign';
    
    %checking if detection far from observation
    reject = [];
    
    for F = 1:numF
        if asign(F) > 0
            reject(F) = est_dist(F,asign(F)) < 50;
        else
            reject(F) = 0;
        end
    end
    
    asign = asign.*reject;
    
    
    %Assign updated detections
    k = 1;
    for F = 1:length(asign)
        if asign(F) > 0
            Q_est(:,k) = Q_est(:,k) + K * (Q_loc_measure(asign(F),:)' - C * Q_est(:,k));
        end
        k = k + 1;
    end
    
    %Update covariance
    P = (eye(4) - K * C) * P;

    %Store Data
    Q_loc_estimateX(t,1:numF) = Q_est(1,1:numF);
    Q_loc_estimateY(t,1:numF) = Q_est(2,1:numF);

    %img = imread(f_list(t).name);
    %imshow(img);
    hold on;
    
    %plot(Y{t}(:), X{t}(:), 'or'); %actual track plot
    
    c_list = ['r' 'b' 'g' 'c' 'm' 'y'];
    set(gca, 'Ydir', 'reverse')
    for Dc = 1:numF
        if ~isnan(Q_loc_estimateX(t,Dc))
            Cz = mod(Dc,6)+1;
            plot(Q_loc_estimateY(t,Dc), Q_loc_estimateX(t,Dc),'o','color',c_list(Cz))
        end
    end
end
    
save('position_estimates.mat','Q_loc_estimateX','Q_loc_estimateY')  


%% Tracking of Dots Returning Radii/Centers 


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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:                                                                        %
% Student Dave: https://www.youtube.com/channel/UCUxiT_SKEUs1oWT6i9P3vPQ             % 
% http://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m-extrema2-m     %
%                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  