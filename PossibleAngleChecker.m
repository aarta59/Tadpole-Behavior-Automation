%just reading and saving movie for later checking 

tadmov = zeros(vidH,vidW,numFrames);
for i = 1:numFrames
    img = read(mov,i);
    img = rgb2gray(img);
    tadmov(:,:,i) = img;
end

[frme, tads] = size(Q_loc_estimateX);

%difference between 1 and 2
X_diff = diff(Q_loc_estimateX);
Y_diff = diff(Q_loc_estimateY);

%difference betweent third frame and second frames
xdiff2 = Q_loc_estimateX(3:length(Q_loc_estimateX),:) - Q_loc_estimateX(2:length(Q_loc_estimateX)-1,:);
ydiff2 = Q_loc_estimateY(3:length(Q_loc_estimateY),:) - Q_loc_estimateY(2:length(Q_loc_estimateY)-1,:);

%angles between second and third frames
tad_angles2 = atan2(xdiff2,ydiff2);
tad_angles2 = abs(tad_angles2);
tad_angles2 = rad2deg(tad_angles2);

logicaltad2 = (tad_angles2 > 1  & tad_angles2 <= 15) | (tad_angles2 >= 165 & tad_angles2 <= 180);
%pads bottom of logical with ones to make it the same size as logical_tad
%assumes that if tad was within angle in 2 frames it is in 3
logicaltad2 = [logicaltad2; ones(1,tads)];

%gets angles between X(2:n,:) - X(1:n-1,:) and Y_diff
tad_angles = atan2(X_diff, Y_diff);



%take absoulte value of angles because direction of movment doesnt matter
%only looking for tadpole 180+/-15 and 0+/-15 (0+/-15 same as abs(15))
tad_angles = abs(tad_angles);

%convert to degrees for simplicity
tad_angles = rad2deg(tad_angles);

%logical array of angles betweent 180+/-15 and 0+/-15 degrees
logical_tad = (tad_angles > 1  & tad_angles <= 15) | (tad_angles >= 165 & tad_angles <= 180);

%within 1 and 2 frames is angle correct?
within2fr = logicaltad2.*logical_tads;


%really just need to know if tadpole is moving within specified angles
%for at least 3 frames might need to interpolate between to get more
%accurate picture of if tadpole is moving correct direction

%loop just shows plot of points along with the actual image of 1 tadpole

for i = 900:numFrames
    imagesc(tadmov(:,:,i+14))

    hold on
    plot(Q_loc_estimateY(i,5),Q_loc_estimateX(i,5),'or')
    plot(Q_loc_estimateY(i+1,5),Q_loc_estimateX(i+1,5),'og')
    plot(Q_loc_estimateY(i+2,5),Q_loc_estimateX(i+2,5),'oc')
    title(['frame' num2str(i+14)])
    set(gca,'YDir','reverse')
   % axis([1 1344 1 1024])
    pause
end


