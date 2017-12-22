% im = read(mov,190);
% im = rgb2gray(im);
% imshow(im)
%hold on 
%plot(Y{190},X{190},'og')
for i = 90:length(Q_loc_estimateX)
    
dz = 4*2;
zx = Q_loc_estimateY(i) - 4;
zy = Q_loc_estimateX(i) - 4;

mn = rectangle('Position',[zx zy dz dz],...
    'Curvature',[1,1],'FaceColor',[0,0,0]);
set(gca,'Ydir','reverse')
axis([1 1344 1 1024])
axis off
imgtad = getframe(gca);
imgtad = frame2im(imgtad);
imgtad = rgb2gray(imgtad);
fulltad(:,:,i) = imgtad;
clf
end

%%

for i = 90:length(Q_loc_estimateX)
%     mov_img = read(mov,i+14);
%     mov_img = rgb2gray(mov_img);
%     imshow(mov_img)
    hold on
    %plot(Q_loc_estimateY(i),Q_loc_estimateX(i),'og')
    d = allradius{i+14}*2;
    px = allcenter{i+14}(:,1) - allradius{i+14};
    py = allcenter{i+14}(:,2) - allradius{i+14};
    for j = 1:length(d)
        h = rectangle('Position',[px(j) py(j) d(j) d(j)],...
            'Curvature',[1,1],'FaceColor',[0,0,0]);
        
    end
    set(gca,'Ydir','reverse')
    axis([1 1344 1 1024])
    axis off
    imgdot = getframe(gca);
    imgdot = frame2im(imgdot);
    imgdot = rgb2gray(imgdot);
    fullimgdot(:,:,i) = imgdot;
    
    clf
end



%%

% mtchpix = (double(imgdot) - double(imgtad)) == 0;
% imshow(mtchpix)
for i = 90:length(Q_loc_estimateX)
co_relate = corr2(fullimgdot(:,:,i),fulltad(:,:,i));
imshowpair(fullimgdot(:,:,i),fulltad(:,:,i))
if co_relate > 0
    encounter = true
else
    non_encounter = false
end
pause
end