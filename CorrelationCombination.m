[fnumber, tadnumber] = size(Q_loc_estimateX);

set(gcf,'Visible','off')
%Drawing dots from center and radii data
frame = zeros(fnumber-89,1);
encount = false(fnumber-89,tadnumber);
for i = 90:length(Q_loc_estimateX)
%     mov_img = read(mov,i);
%     mov_img = rgb2gray(mov_img);
%     imshow(mov_img)
%     hold on
%     plot(Q_loc_estimateY(i-14,1),Q_loc_estimateX(i-14,1),'og')
    
    
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
    
        co_relate = corr2(imgdot,imgtad);
        %imshowpair(fullimgdot(:,:,i),fulltad(:,:,i-14))
        if co_relate > 0
            encounter = true;
        else
            encounter = false;
        end
    
    
            frame(i-89,1) = i;
            encount(i-89,k) = encounter;
    end
       
       
end

framesAndEncounters = [frame encount];
save('frame_number_and_encounter.mat','frameAndEncounters')
