%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       T.A.D.9000 Visualization                        %
%                        Aaron Ta, October 2019                         %
%              Cline Lab, Dorris Center for Neuroscience                %
%           Scripps Research Institute, La Jolla, California            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for visualizing the T.A.D.9000 function behavior post-analysis.%
%                                                                       %
% REQUIRED FILES IN FOLDER (outputted by T.A.D.9000):                   %
% ~ original video (mp4 by default)                                     %
% ~ position_estimates.mat (estimated tadpole coordinates)              %
% ~ dots_centers_radii.mat (detected dot coordinates)                   %
% ~ encounter_matrix.mat (table of tadpole-dot encounters)              %
% ~ Final_Avoidance_Data.mat (list of valid encounters + turning check) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETECTION KEY:                                                        %
%                                                                       %
% ~ small gray = default/no collision                                   %
% ~ small white = unconsidered collision                                %
% ~ medium white = collision w/ no turn                                 %
% ~ large white = collision w/ turn                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of frames specified in T.A.D.9000 to be discarded as startup
frameStartup = 15;
%Number of frames specified in T.A.D.9000 before dots appear
dotFrameStart = 450;

%FPS of output video (default is 30 fps)
setFrameRate = 30;

%Color/size definitions (color value in grayscale, size in px*2)
white = 255;
gray = 170;
darkTint = 150;

small=2;
medium=4;
large=6;

%small gray = default/no collision
%small white = unconsidered collision
%medium white = collision w/ no turn
%large white = collision w/ turn

%Load original video
mov_list = dir('*mp4'); %Change file extension for .avi
mov_name = mov_list.name;
mov = VideoReader(mov_name);

%Obtain video length, height, and width; initialize a copy
frameLength = uint16(mov.Duration) * uint16(mov.FrameRate);
netFrames = frameLength - frameStartup + 1;
vidHeight = mov.Height;
vidWidth = mov.Width;

tmpVid=int16(zeros(vidHeight,vidWidth,netFrames));
for i = 1:netFrames
    tmpVid(:,:,i) = int16(rgb2gray(read(mov,i+frameStartup-1)));
end


%Load the position estimate file
pos_estimate=load('position_estimates.mat');
%X- and Y- position of each tadpole at each frame
xPos=double(int16(pos_estimate.Q_loc_estimateX));
yPos=double(int16(pos_estimate.Q_loc_estimateY));
xPos=num2cell(transpose(xPos),1);
yPos=num2cell(transpose(yPos),1);

%Load the dot centers and radii file
centers_radii=load('dot_centers_radii.mat');
%Coordinates of dot centers and radii of dots
centers=centers_radii.allcenter;
radii=centers_radii.allradius;

dotFrameStart = 1;
while size(centers{dotFrameStart}) == 0
    dotFrameStart = dotFrameStart + 1;
end

%Convert center/radii information into doubles
for r = 1:length(centers)
    centers{r}=double(int16(centers{r}));
end

for r = 1:length(radii)
    radii{r}=double(int16(radii{r}));
end

%Load the encounter matrix file
enc_matrix=load('encounter_matrix.mat');
%Table of whether or not each tadpole encountered a dot at each frame
enc_matrix=enc_matrix.encounterMatrix;

%Load the avoidance file
avoid_matrix=load('Final_Avoidance_Data.mat');
%Table of whether or not each valid encounter also involved a turn
avoid_matrix=avoid_matrix.dataStore;

clear mov;
clear centers_radii;
clear pos_estimate;

display(newline)
display('All files loaded.');

%Annotate all frames
%Load mesh grid of video for circle creation
[columnsInImage, rowsInImage] = meshgrid(1:vidWidth, 1:vidHeight);
for f = 1:netFrames
    %Approximate dot detections with dark gray squares
    for i = 1:length(centers{f+frameStartup-1})
        tmpRadii = radii{f+frameStartup-1}(i);
        xCenter = centers{f+frameStartup-1}(i,2);
        yCenter = centers{f+frameStartup-1}(i,1);
        circlePixels = (rowsInImage - xCenter).^2 + (columnsInImage - yCenter).^2 <= tmpRadii.^2;
        tmpImg = tmpVid(:,:,f);
        tmpImg(circlePixels) = 50;
        tmpVid(:,:,f) = tmpImg;
%            for p = -radii{f+frameStartup-1}(i):radii{f+frameStartup-1}(i)
%               for q = -radii{f+frameStartup-1}(i):radii{f+frameStartup-1}(i)
%                   xSize = centers{f+frameStartup-1}(i,2) + p;
%                   ySize = centers{f+frameStartup-1}(i,1) + q;
%                   if xSize < 1
%                       xSize = 1;
%                   end
%                   if ySize < 1
%                       ySize = 1;
%                   end
%                   darken = tmpVid(xSize,ySize,f) - darkTint;
%                   if darken < 0
%                       darken = 0;
%                   end
%                   tmpVid(xSize,ySize,f) = darken;
%               end
%            end
    end
    %Approximate tadpole position estimations with gray squares that turn
    %white upon a detected collision
    for i = 1:length(xPos{f})
        if f > (dotFrameStart-frameStartup) && enc_matrix(f + frameStartup - dotFrameStart,i) == 1
            %Make a white box on estimated position on frames w/ collisions
            color=white;
        else
            %Make a gray box on estimated position otherwise
            color=gray;
        end
        if ismember(f + frameStartup - 1,avoid_matrix{1,i}(:,1))
            if avoid_matrix{1,i}(find(avoid_matrix{1,i}(:,1)==(f + frameStartup - 1)),2) == 1
                %Large box for valid encounter + turning
                dsize=large;
            else
                %Small box for valid encounter w/o turning
                dsize=small;
            end
        else
            %Medium box for unconsidered encounter
            dsize=medium;
        end
        xCenter = xPos{f}(i);
        yCenter = yPos{f}(i);
        circlePixels = (rowsInImage - xCenter).^2 + (columnsInImage - yCenter).^2 <= dsize.^2;
        tmpImg = tmpVid(:,:,f);
        tmpImg(circlePixels) = color;
        tmpVid(:,:,f) = tmpImg;
%        for p = -dsize:dsize
%            for q = -dsize:dsize
%               xSize = xPos{f}(i) + p;
%               ySize = yPos{f}(i) + q;
%               if xSize < 1
%                   xSize = 1;
%               elseif xSize > vidHeight
%                   xSize=vidHeight;
%               end
%               if ySize < 1
%                   ySize = 1;
%               elseif ySize > vidWidth
%                   ySize=vidWidth;
%               end
%               tmpVid(xSize,ySize,f) = color;
%            end
%        end
    end
    %Display checkpoint message once 50% of frames are processed
    if (mod(f,2) == 0 & (f == f/2)) | (mod(f,2) == 1 & ((f-1) == netFrames/2))
      display('50% of frames processed.');
    end
end

display('100% of frames processed.');
display('Dot detections approximated.');
display('Tadpole collisions annotated.');

%Write modified video to file        
outVid2 = VideoWriter(fullfile(strcat(extractBefore(mov_name,'.'),'_AnnotatedCircleSlow.mp4')));
outVid2.FrameRate=setFrameRate;
open(outVid2);
for f = 1:netFrames
    writeVideo(outVid2,double(tmpVid(:,:,f))/256)
    if (mod(f,2) == 0 & (f == f/2)) | (mod(f,2) == 1 & ((f-1) == netFrames/2))
      display('50% of frames written to file.');
    end
end
close(outVid2)


display('100% video output written.');
display('Done.');
display(newline)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code References:                                                                   %
% MATLAB FAQ: https://matlab.fandom.com/wiki/FAQ                                     %                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
