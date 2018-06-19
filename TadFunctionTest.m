%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             T.A.D.9000                                %
%                 Christopher Marshall, November 2017                   %
%              Cline Lab, Dorris Center for Neuroscience                %
%           Scripps Research Institute, La Jolla, California            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [encAvg,numEncount,numAvoid] = TadFunctionTest(mov)
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
sigmah = 6;   
h = fspecial('log', hsizeh, sigmah);

%% Iteratively finding tadpoles from blobs

%enter number of tadpoles in video for auto threshold
%tad_number = 6;

%Starting frame for detection (start at 15 for even brightness)
s_frame = 15;

%Initialize cells for detection coordinates
X = cell(1,numFrames-(s_frame-1)); %detection X coordinate 
Y = cell(1,numFrames-(s_frame-1));  %detection Y coordinate 

for i = 1:numFrames-(s_frame-1)
    %img_real = (read(mov,i));
    bck_img = double(bck_img);
    img = noDot_img(:,:,i+(s_frame-1)); 
    sub_img = (img - bck_img);

    %Blob Filtering
    blob_img = conv2(sub_img,h,'same');
 
    %Thresholding level for blob (0.032 ex 95, 0.4 ex 96)
    idx = find(blob_img < 0.1); 
    blob_img(idx) = nan;
    
    %Finds peak indices for blobs
    [~,imax,~,~] = extrema2(blob_img); 
    

    %Section for auto threshold based on number of expected detections
%     thresh = 0.040;
%     while length(imax) > tad_number 
%         idx = find(blob_img < thresh); 
%         blob_img(idx) = nan;
%         [zmax,imax,zmin,imin] = extrema2(blob_img);
%             
%         if length(imax) > tad_number
%             thresh = thresh + 0.002;
%         end
%     end
    
    [X{i},Y{i}] = ind2sub(size(blob_img),imax);
    
    if i == (round((numFrames-(s_frame-1))/2))
        disp('Tadpole detection is 50% complete')
    elseif i == (numFrames-(s_frame-1))
        disp('Tadpole detection is 100% complete')  
    end
    
    %Plot of raw detections with threshold overlay
%      imagesc(blob_img)
%     hold on
%     for j = 1:length(X{i})
%        plot(Y{i}(j),X{i}(j),'or')
%     end
%     axis off
%     
%      pause
%      
end

save('raw_tad_detections.mat','X','Y')

%% Kalman Filter Variable Definitions 

%These values should be changed to allow for more accurate assignment
%Current values are working: (dt=1,u=0,tnm=1,tmnx=0.5,tnmy=0.5)
dt = 1; %sampling rate
u = 0; %starting acceleration magnitude 
Tad_noise_mag = 3; %variability in tadpole speed 
tmn_x = 2; %noise in horizontal direction, x-axis
tmn_y = 2; %noise in vertical direction, y-axis

%Process noise into covariance matrix (Ex)
Ez = [tmn_x 0; 0 tmn_y];
Ex = [dt^4/4 0 dt^3/2 0; 0 dt^4/4 0 dt^3/2; dt^3/2 0 dt^2 0; 0 dt^3/2 0 dt^2].*Tad_noise_mag^2;

%Estimate of initial position variance (covariance matrix)
P = Ex;

%2-D updates for coefficent matrix
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
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
    
    [asign, ~] = munkres(est_dist);
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

%% Removing Bad Detections

%removes columns where detection is not in frame
removeX = any(Q_loc_estimateX<0 | Q_loc_estimateX>vidH);
removeY = any(Q_loc_estimateY<0 | Q_loc_estimateY>vidW);
removePoints = [removeX; removeY];
removePoints = any(removePoints == 1);

Q_loc_estimateX(:,removePoints) = [];
Q_loc_estimateY(:,removePoints) = [];

%% Tracking of Dots Returning Radii/Centers 

%Cropping movie to remove false dot recognition
% ytop = 140;
% ybott = 910;
% xleft = 90;
% xright = 1250;

%try with these new values
% ytop = 160;
% ybott = 900;
% xleft = 100;
% xright = 1230;

%values for channel system (xright = 1280)
ytop = 130;
ybott = 900;
xleft = 80;
xright = 1125;

% prealocate cells for dot position storage
allcenter = cell(1,numFrames); %X,Y coordinate centers of dots index
allradius = cell(1,numFrames); %radius of each detected dot 
dotzeros = uint8(zeros(vidH,vidW));
dotLoopLength = length(90:numFrames);

for i = 90:numFrames
    orig_img = read(mov,i);
    orig_img = rgb2gray(orig_img);
    croped_orig = orig_img(ytop:ybott,xleft:xright,:);
    %backg = imopen(croped_orig, strel('disk',23));
    backg = uint8(bck_img(ytop:ybott,xleft:xright,:));
    minus_bck = croped_orig - backg;
    adj_mius = imadjust(minus_bck);
    str_mius = imopen(adj_mius, strel('disk',4));
    %str_mius = str_mius*2;
    dotzeros(ytop:ybott,xleft:xright) = str_mius;

    [allcenter{i}, allradius{i}] = imfindcircles(dotzeros,[10 20],...
        'ObjectPolarity','bright', 'Sensitivity',0.90);
    
    if i == (round(dotLoopLength*0.25))
        disp('Dot detection is 25% complete')
    elseif i == (round(dotLoopLength*0.50))
        disp('Dot detection is 50% complete')
    elseif i == (round(dotLoopLength*0.75))
        disp('Dot detection is 75% complete')
    elseif i == numFrames
        disp('Dot detection is 100% complete')
    end
    
end

%Save radii and centers of dots
save('dot_centers_radii.mat','allcenter','allradius')

%% Drawing dots and tadpoles then computing correlation 
%Clipping X and Y positions so at index 1, frame is 90
clippedX = Q_loc_estimateX(76:end,:);
clippedY = Q_loc_estimateY(76:end,:);

%Clipping dot centers and radius so index 1, frame is 90
clipCenters = allcenter(90:end);
clipRadius = allradius(90:end);

[frme, tads] = size(clippedX);

%difference tadpole is currently at 
xdiff1 = clippedX(1:length(clippedX)-1,:) - clippedX(2:length(clippedX),:);
ydiff1 = clippedY(1:length(clippedY)-1,:) - clippedY(2:length(clippedY),:);

%gets angles between X(2:n,:) - X(1:n-1,:) and Y_diff
tad_angles1 = atan2d(xdiff1, ydiff1);

%take absoulte value of angles because direction of movment doesnt matter
%only looking for tadpole 180+/-15 and 0+/-15 (0+/-15 same as abs(15))
tad_angles1 = abs(tad_angles1);

xdiff2 = diff(clippedX);
ydiff2 = diff(clippedY);

%angles between second and third frames
tad_angles2 = atan2d(xdiff2, ydiff2);
tad_angles2 = abs(tad_angles2);

%padding last value of tadpole angles with 0
t_pad_angle = [tad_angles2; zeros(1,tads)];

%points around tadpole to plot
theta = (0:360)';

%preallocate encounter matrix
encounterMatrix = zeros(frme,tads);

for i = 1:frme
    
    cenXYRad = [clipCenters{i} clipRadius{i}];
    
    %loop starts with column 1 == tadpole 1
    for j = 1:tads
        
        %guess of r=15 pixels (eyes 13 pixels from gut of tadpole)
        pline_x = 15*cos(theta) + clippedX(i,j);
        pline_y = 15*sin(theta) + clippedY(i,j);
        
        %looks only at points in half-circle in direction of movment
        if t_pad_angle(i,j) < 90
            tad_visionIdx = find(pline_y > clippedY(i,j));
            pline_x = pline_x(tad_visionIdx);
            pline_y = pline_y(tad_visionIdx);
        else
            tad_visionIdx = find(pline_y < clippedY(i,j));
            pline_x = pline_x(tad_visionIdx);
            pline_y = pline_y(tad_visionIdx);
        end
        
        %looks at every dot center point and finds distance from dot center
        %to all points around tadpole
        dotEncount = 0;
        for k = 1:length(cenXYRad)

            dist_Circ = sqrt((pline_x - cenXYRad(k,2)).^2 + (pline_y - cenXYRad(k,1)).^2);
            
            if any(dist_Circ <= cenXYRad(k,3))
                dotEncount = dotEncount + 1;
            end
        end
        
        if dotEncount > 0
            encounterMatrix(i,j) = 1;
        else
            encounterMatrix(i,j) = 0;
        end
        
    end
end

save('encounter_matrix.mat','encounterMatrix')


%% Logic for angle checking and velocity checking

%section does two checks for tadpole movment. first check is for the angle
%of the tadpole between frame "n" and "n-1" the second check is for the
%angle between "n" and "n+1" to see if eye is perpendicular to dot

%logical array of angles between 180+/-15 & 0+/-15 degrees for frame (1-2)
logicaltad1 = ((tad_angles1 >= 0)  & (tad_angles1 <= 15)) | ((tad_angles1 >= 165) & (tad_angles1 <= 180));

%logical array of angles between 180+/-15 & 0+/-15 degrees for frame (2-3)
logicaltad2 = ((tad_angles2 >= 0) & (tad_angles2 <= 15)) | ((tad_angles2 >= 165) & (tad_angles2 <= 180));

%within 1 and 2 frames is angle correct?
within2fr = logicaltad2.*logicaltad1;

%checks velocity of tadpole
f_rate = 15; %f/s
t_disp = 1/f_rate; %sec (time between frames)

Vx = xdiff1./t_disp;
Vy = ydiff1./t_disp;

%velocity [=] pixels/second
Vtot = sqrt(Vx.^2 + Vy.^2);

%logic check for velocity less than 10000 and greater than 50 pix/sec
%50 pix/s seems like the threshold for real movment

velLogic = (Vtot > 50) & (Vtot < 10000);

%here is if the angle for 3 frames is ok and the velocity is within range
within2frAndVelocity = velLogic.*within2fr;

%pads bottom with zeros assumes last frame velocity and angles are 0
%within2frAndVelocity = [within2frAndVelocity; zeros(1,tads)];
within2frAndVelocity = [zeros(1,tads); within2frAndVelocity];

actualEncounters = encounterMatrix.*within2frAndVelocity;

actualFramesAndEncount = [((1:frme)+89)' actualEncounters];

%removes last 8 encounters
remFrames = any(actualFramesAndEncount((frme-8):frme(end),2:end) == 1);
remIdx = find(remFrames == 1);
actualFramesAndEncount((frme-8):frme(end),(remIdx+1)) = 0;

%% Event Detection 

%need to know if tadpoles angle changed between 90 and 180 degrees of where
%encounter occured within 8 frames
% look at frame where encounter happened + 8 frames 

for i = 2:length(actualFramesAndEncount(1,:))
    
    encounter = find(actualFramesAndEncount(:,i) == 1);
    
    %point where encounter occurs
    enc_pointX = clippedX(encounter,i-1);
    enc_pointY = clippedY(encounter,i-1);

    %points from encounter frame (n) to 8 frames in future
    future_points = encounter+(1:8);
    [r,c] = size(future_points);
    
    fut_pointX = zeros(r,c);
    fut_pointY = zeros(r,c);
    
    for k = 1:c
        fut_pointX(:,k) = clippedX(future_points(:,k),i-1);
        fut_pointY(:,k) = clippedY(future_points(:,k),i-1);
    end

    %both encounter point and future points
    a = [enc_pointX fut_pointX]';
    b = [enc_pointY fut_pointY]';

    %difference between encounter point (n) and future point (n+1)
    event_diffX = diff(a);
    event_diffY = diff(b);

    %angle between encounter frame (n) and future frame (n+(1:8))
    event_angle = atan2d(event_diffX,event_diffY);
    event_angle = abs(event_angle);

    [nrow, ncol] = size(event_angle);
    event_count = zeros(nrow,ncol);
    
    for c = 1:ncol
        for r = 1:nrow
            
            %increase tolerance for event here to do this
            %lower first if statment event_angle(r,c) > 85-tolerance
            %increase elseif statment event_angle(r,c) <= 95+tolerance
            if event_angle(1,c) < 90 && event_angle(r,c) >= 70
                event = true;
            elseif event_angle(1,c) > 90 && (event_angle(r,c) <= 110 && event_angle(r,c) >= 0)
                event = true;
            else
                event = false;
            end
            event_count(r,c) = event;
        end
    end
    
    event = any(event_count)';
    encounterFrameandEvent = [(encounter+89) event];
    
    dataStore{i-1} = encounterFrameandEvent;
end

save('Final_Avoidance_Data.mat','dataStore')
%END OF PROGRAM

%VALUE of 'dataStore' contains frame at which encounter occurred and if that
%encounter triggered an event 

%% Provide Average 
encAvg = zeros(1,length(dataStore));
numEncount = zeros(1,length(dataStore));
numAvoid = zeros(1,length(dataStore));

for i = 1:length(dataStore)
    sumOnes = sum(dataStore{i}(:,2) == 1);
    avg = sumOnes/length(dataStore{i});
    encAvg(i) = avg;
    numEncount(i) = length(dataStore{i});
    numAvoid(i) = sumOnes;
end

tab = table(encAvg', numAvoid', numEncount', 'VariableNames', {'AvoidanceIndex',...
    'NumberAvoidances', 'NumberEncounters'});
writetable(tab)
end

%% References and Dependencies 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code References:                                                                   %
% Student Dave: http://studentdavestutorials.weebly.com/                             % 
%                                                                                    %        
%                                                                                    %
% Runtime Dependencies:                                                              %
% https://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m--extrema2-m   %                                                                                %
% https://www.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker          %                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  