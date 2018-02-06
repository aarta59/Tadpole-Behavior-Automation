%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             T.A.D.9000                                %
%                 Christopher Marshall, November 2017                   %
%              Cline Lab, Dorris Center for Neuroscience                %
%           Scripps Research Institute, La Jolla, California            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Choose location of video and file to be run
directory = uigetdir;
cd(directory);
moviename = uigetfile('*.mov');
folder = fullfile(directory);
movFullFile = fullfile(folder, moviename);
mov = VideoReader(movFullFile);

%Video dimentions
numFrames = mov.NumberOfFrames;
vidH = mov.Height;
vidW = mov.Width;

%Initializations 
img = zeros(vidH,vidW,numFrames);
noDot_img = zeros(vidH,vidW,numFrames);

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

%% Iteratively finding tadpoles from blobs

%Starting frame for detection (start at 15 for even brightness)
s_frame = 15;

%Initialize cells for detection coordinates
X = cell(1,numFrames-(s_frame-1)); %detection X coordinate 
Y = cell(1,numFrames-(s_frame-1));  %detection Y coordinate 

for i = 1:numFrames-(s_frame-1)-500
    %img_real = (read(mov,i));
    bck_img = double(bck_img);
    img = noDot_img(:,:,i+(s_frame-1)); 
    sub_img = (img - bck_img);

   
    
    %Blob Filtering
    blob_img = conv2(sub_img,h,'same');
 
    %Thresholding level for blob
    idx = find(blob_img < 0.032); 
    blob_img(idx) = nan;
    
   
    %Finds peak indices for blobs
    %[zmax,imax] = max(blob_img(:)); %ONLY WORKS FOR 1 TADPOLE CURRENTLY
    [zmax,imax,zmin,imin] = extrema2(blob_img); %WORKS FOR MULTIPLE TADS
    
    [X{i},Y{i}] = ind2sub(size(blob_img),imax);
    
    %Plot of raw detections with threshold overlay
%     imagesc(blob_img)
%     hold on
%     for j = 1:length(X{i})
%        plot(Y{i}(j),X{i}(j),'or')
%     end
%     axis off
%     
%     pause
     
end

save('raw_tad_detections.mat','X','Y')

%% Kalman Filter Variable Definitions 

%These values should be changed to allow for more accurate assignment
%Current values are working: (dt=1,u=0,tnm=2,tmnx=0.05,tnmy=0.05)
dt = 1; %sampling rate
u = 0; %starting acceleration magnitude 
Tad_noise_mag = 2; %variability in tadpole speed 
tmn_x = 0.05; %noise in horizontal direction, x-axis
tmn_y = 0.05; %noise in vertical direction, y-axis

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
Q = [X{1} Y{1} zeros(length(X{1}),1) zeros(length(X{1}),1)]';
Q_est = nan(4,2000);
Q_est(:,1:size(Q,2)) = Q; %initial location estimate
Q_loc_estimateY = nan(2000);
Q_loc_estimateX = nan(2000);
numDet = size(X{1},1); %number of detections
numF = find(isnan(Q_est(1,:))==1, 1)-1; %number of estimates

for t = 1:length(X)
    
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
    K = (P * C')/(C * P * C' + Ez);
    
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
            reject(F) = est_dist(F,asign(F)) < 100;
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
    
    %Remove nan values and store 
    Q_loc_estimateX(isnan(Q_loc_estimateX)) = [];
    Q_loc_estimateY(isnan(Q_loc_estimateY)) = [];
     

%     imagesc(noDot_img(:,:,t));
%     hold on;
%     
%     %plot(Y{t}(:), X{t}(:), 'or'); %actual track plot
%     
%     c_list = ['r' 'b' 'g' 'c' 'm' 'y'];
%     set(gca, 'Ydir', 'reverse')
%     for Dc = 1:numF
%         if ~isnan(Q_loc_estimateX(t,Dc))
%             Cz = mod(Dc,6)+1;
%             plot(Q_loc_estimateY(t,Dc), Q_loc_estimateX(t,Dc),'o','color',c_list(Cz))
%         end
%     end
%    
end

save('position_estimates.mat','Q_loc_estimateX','Q_loc_estimateY')  

%% Plots location estimates

%get dimentions of location estimates
[numPositions, numDetections] = size(Q_loc_estimateX);

figure
hold on
axis off
set(gca,'YDir','reverse')
c_list = ['r' 'b' 'g' 'c' 'm' 'y'];
for j = 77:numPositions-500
    for i = 1:numDetections
        cz = mod(i,6)+1;
        plot(Q_loc_estimateY(j,i),Q_loc_estimateX(j,i), 'o', 'color', c_list(cz))
    end
    pause
end
    

%% Tracking of Dots Returning Radii/Centers 

%Cropping movie to remove false dot recognition
ytop = 140;
ybott = 910;
xleft = 90;
xright = 1250;

% prealocate cells for dot position storage
allcenter = cell(1,numFrames); %X,Y coordinate centers of dots index
allradius = cell(1,numFrames); %radius of each detected dot 
dotzeros = uint8(zeros(vidH,vidW));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only reason to reoverlay image on zeros is to show %
% the dots on the original image. In final revision  %
% remove the dotzeros part and move straight to next %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 90:numFrames
    orig_img = read(mov,i);
    orig_img = rgb2gray(orig_img);
    croped_orig = orig_img(ytop:ybott,xleft:xright,:);
    backg = imopen(croped_orig, strel('disk',23));
    minus_bck = croped_orig - backg;
    adj_mius = imadjust(minus_bck);
    str_mius = imopen(adj_mius, strel('disk',4));
    str_mius = str_mius*2;
    dotzeros(ytop:ybott,xleft:xright) = str_mius;

    [allcenter{i}, allradius{i}] = imfindcircles(dotzeros,[10 20],...
        'ObjectPolarity','bright', 'Sensitivity',0.90);    
end

%Save radii and centers of dots
save('dot_centers_radii.mat','allcenter','allradius')


%% Drawing dots and tadpoles then computing correlation 

%sizing location matrix (index starts at frame 15)
[fnumber, tadnumber] = size(Q_loc_estimateX);

%plots not visible
set(gcf,'Visible','off')

%Drawing dots from center and radii data
% frame = zeros(numFrames-89,1);
% encount = false(numFrames-89,tadnumber);

for i = 90:numFrames
    mov_img = read(mov,i);
    mov_img = rgb2gray(mov_img);
    imshow(mov_img)
    hold on
    plot(Q_loc_estimateY(i-14,2),Q_loc_estimateX(i-14,2),'og')

    d = allradius{i}*2;
    px = allcenter{i}(:,1) - allradius{i};
    py = allcenter{i}(:,2) - allradius{i};
    for j = 1:length(d)
        h = rectangle('Position',[px(j) py(j) d(j) d(j)],...
            'Curvature',[1,1],'FaceColor',[0,0,0]);  
    end
    set(gca,'Ydir','reverse')
    axis([1 1344 1 1024])
    axis off
    imgdot = getframe(gcf);
    imgdot = frame2im(imgdot);
    imgdot = rgb2gray(imgdot);
    clf
    
    for k = 1:tadnumber
        dz = 4*2;
        zx = Q_loc_estimateY(i-14,k) - 4;
        zy = Q_loc_estimateX(i-14,k) - 4;

        mn = rectangle('Position',[zx zy dz dz],...
            'Curvature',[1,1],'FaceColor',[0,0,0]);
        set(gca,'Ydir','reverse')
        axis([1 1344 1 1024])
        axis off
        imgtad = getframe(gcf);
        imgtad = frame2im(imgtad);
        imgtad = rgb2gray(imgtad);
        clf
        
        %corr2 shows when tadpole crosses dot
        co_relate = corr2(imgdot,imgtad);
        
        % negative correlation means dot and tadpole are far
        % positive correlation means dot and tadpole intersect
        if co_relate > 0
            encounter = true;
        else
            encounter = false;
        end

            frame(i-89,1) = i;
            encount(i-89,k) = encounter;

    end

end

%index for frames and encounters starts at frame 90 (where dots begin)

framesAndEncounters = [frame encount];
save('frame_number_and_encounter.mat','frameAndEncounters')

%%%%%%%%%%%%%%%%%-NOTES-%%%%%%%%%%%%%%%%%
% location estimates are frame-14 off   %
% of dot drawing due to re-indexing of  %
% locations starting at 1 (prev at 15)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section could be used to find location of intersection %
%                                                             %
% mtchpix = (double(imgdot) - double(imgtad)) == 0;           % 
% imshow(mtchpix)                                             %
% clear frame encount                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Visual Check if encounter occured

% Below is just a simple check to see if the encounter actually did occur
for i = 90:length(Q_loc_estimateX)
    mov_img = read(mov,i);
    mov_img = rgb2gray(mov_img);
    imshow(mov_img)
    hold on
    plot(Q_loc_estimateY(i-14,1),Q_loc_estimateX(i-14,1),'og')
    title(['frame' num2str(i)])
    pause
    
end

%% Logic for angle checking and velocity checking

%just reading and saving movie for later checking 

tadmov = zeros(vidH,vidW,numFrames);
for i = 1:numFrames
    img = read(mov,i);
    img = rgb2gray(img);
    tadmov(:,:,i) = img;
end
tadmov = uint8(tadmov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clippedX = Q_loc_estimateX(76:end,:);
clippedY = Q_loc_estimateY(76:end,:);

[frme, tads] = size(clippedX);

%difference between 1 and 2 frames
xdiff1 = diff(clippedX);
ydiff1 = diff(clippedY);

%gets angles between X(2:n,:) - X(1:n-1,:) and Y_diff
tad_angles1 = atan2d(xdiff1, ydiff1);

%take absoulte value of angles because direction of movment doesnt matter
%only looking for tadpole 180+/-15 and 0+/-15 (0+/-15 same as abs(15))
tad_angles1 = abs(tad_angles1);


%logical array of angles betweent 180+/-15 and 0+/-15 degrees
logicaltad1 = ((tad_angles1 > 1)  & (tad_angles1 <= 15)) | ((tad_angles1 >= 165) & (tad_angles1 <= 195));

%difference between 2 and 3 frames
xdiff2 = clippedX(3:length(clippedX),:) - clippedX(2:length(clippedX)-1,:);
ydiff2 = clippedY(3:length(clippedY),:) - clippedY(2:length(clippedY)-1,:);

%angles between second and third frames
tad_angles2 = atan2d(xdiff2, ydiff2);
tad_angles2 = abs(tad_angles2);


logicaltad2 = ((tad_angles2 > 1) & (tad_angles2 <= 15)) | ((tad_angles2 >= 165) & (tad_angles2 <= 195));

%pads bottom of logical with ones to make it the same size as logical_tad
%assumes that if tad was within angle in 2 frames it is in 3
logicaltad2 = [ones(1,tads); logicaltad2];

%within 1 and 2 frames is angle correct?
within2fr = logicaltad2.*logicaltad1;

%checks velocity of tadpole
f_rate = 15; %f/s
t_disp = 1/f_rate; %sec (time between frames)


Vx = xdiff1./t_disp;
Vy = ydiff1./t_disp;

%velocity [=] pixels/second
Vtot = sqrt(Vx.^2 + Vy.^2);

%logic check for velocity less than 10000 and greater than 10 pix/sec
velLogic = (Vtot > 10) & (Vtot < 10000);

%here is if the angle for 3 frames is ok and the velocity is within range
within2frAndVelocity = velLogic.*within2fr;

%pads bottom with zeros assumes last frame velocity and angles are 0
within2frAndVelocity = [within2frAndVelocity; zeros(1,tads)];

actualEncounters = frameAndEncounters(:,2:end).*within2frAndVelocity;

actualFramesAndEncount = [frameAndEncounters(:,1) actualEncounters];

%really just need to know if tadpole is moving within specified angles
%for at least 3 frames might need to interpolate between to get more
%accurate picture of if tadpole is moving correct direction

%loop just shows plot of points along with the actual image of 1 tadpole
whichone = 2;


for i = 1:frme-500
    imshow(tadmov(:,:,i+89))
    hold on
        for k = 1:tads
            plot(clippedY(i,k),clippedX(i,k),'or')
            plot(clippedY(i+1,k),clippedX(i+1,k),'og')
            plot(clippedY(i+2,k),clippedX(i+2,k),'oc')
            %title(['frame: ' num2str(i+103) ' encounter?: ' num2str(actualEncounters(i,whichone)) ' angle?: ' num2str(within2frAndVelocity(i,whichone))])
        end
    set(gca,'YDir','reverse')
    axis off
   % axis([1 1344 1 1024])
    pause
    
end



%% TO DO LIST
%               MOST IMPORTANT
% 1)check for angle of tadpole parallel to X (+/- 15 deg)(DONE)

%               Later Work
% 2) compute velocity of tadpole over range of video (DONE)
% 3) way to correlate matrices without re-saving images from figures
% 4) save image and frame number of encounter for later checking (DONE)
% 5) project average position for tadpole eyes instead of tracking gut
%   5a) dont need this if you know direction tadpole moves in


%% Plot circles instead of drawing

% [This will most likely be deleted in final revision] 

% functional, could be used instead of correlation method
% by finding if tadpole comes within certain distance from
% the computed points below
% possibly faster method will try to get working after

% theta = 0:0.01:(2*pi);
% [pline_x] = allradius{i}*cos(theta) + allcenter{i}(:,1);
% [pline_y] = allradius{i}*sin(theta) + allcenter{i}(:,2);
% plot(pline_x,pline_y,'.r');
% set(gca,'Ydir','reverse');
% axis off







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:                                                                        %
% Student Dave: http://studentdavestutorials.weebly.com/                             % 
%                                                                                    %
%                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  