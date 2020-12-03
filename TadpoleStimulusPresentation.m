clear;
%Suppress PsychToolbox Warnings
Screen('Preference', 'SkipSyncTests', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% set CCD prameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switchinterval = 0; %wait seconds between switch drifting side
waitsec = 10; % waiting time before showing stimulus, default:10 s
vscreen = 18.45; % verticle tank dimension (actual size of projection on podium in cm)
hscreen = 13.84; % horizontal tank dimension (actual size of projection on podium in cm)
vres = 800; % verticle screen resolution
hres = 600; % horizontal screen resolution
nrowOffset = 0; %REVERT!
scaleFactor = 0; %TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up camera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vidobj = videoinput('hamamatsu',1);
%28.1fps capture rate x (60s + 20s startup) ~= 2250 frames
vidobj.FramesPerTrigger = 2250;
vidobj.LoggingMode = 'disk';
src = getselectedsource(vidobj);
%src.BacklightCompensation = 'on';
%src.ExposureMode = 'auto';
set(src, 'ExposureTime', 0.006)
%src.WhiteBalanceMode = 'auto';
%src.Exposure = -6;
%src.Brightness = -60;
%src.Gamma = 165;
%src.Contrast = 0;
%src.Hue=0;
%src.Saturation=0;
%src.WhiteBalance=4600;

% Select a codec for the AVI file.
%logfile.Compression = 'MSVC';
% Since grayscale images will be acquired, a colormap is required.
%logfile.Colormap = gray(256);
vidobj.ReturnedColorSpace = 'grayscale';
% Set region of interest to a 800x600 field with an offset
vidobj.ROIPosition = [120 248 800 600];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select file name and path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Pick a Directory to Save Your Image and Data'
prompt={'Enter the file path to save your image file from CCD and parameters)',...
    'Enter the file name:'};
name='Input the file name and path for saving files';
numlines=1;
defaultanswer={'C:\Users\Cline Lab\Desktop\20201021_AT',...
    'test-'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
filefolder=fullfile(answer{1},answer{2});

if (exist(filefolder) == 0)
    mkdir(filefolder);
end

cd(filefolder)

vidname=char(strcat(strcat(strcat(answer(1),'\'),answer(2)),'.mp4'))
% Create an AVI file object.
logfile = VideoWriter(vidname,'MPEG-4')
% Set frame rate of AVI file
%logfile.FrameRate=15;
logfile.Quality=95;
% Configure the video input object to use the AVI file object.
vidobj.DiskLogger = logfile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select type of stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = { 'Random dots','Full Gratings','Single Gratings'};
[sel,ok] = listdlg('PromptString','Select a mode:',...
    'SelectionMode','single',...
    'ListSize',[150 100],...
    'Name', 'Mode',...
    'ListString',str);
if ok==0
    error('No mode selected!')
end
str2 = {'Half image', 'Entire image'};
[sel2,ok2] = listdlg('PromptString','Select a mode:',...
    'SelectionMode','single',...
    'ListSize',[150 50],...
    'Name', 'Mode2',...
    'ListString',str2);
if ok2==0
    error('No mode selected!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotorgrid={'dots' 'grid' 'grid'};
%Default Parameters for screen size,resolution,velocity,refresh,size
di1={'2','16','45','100','0.2','49','30'};
% the number of dots is defined by 1.4 dots/cm2 from Carlos paper
prompt  = {'Velocity (cm/s)','Stimulus refresh rate (ms)',...
    'Duration of each cycle (s)', 'Number of dots', ...
    'Circle radius (cm),Rang(0-0.6)',...
    'Xenopus stage','Bar duration',};
defAns={char(di1(1)),char(di1(2)),...
    char(di1(3)), char(di1(4)),...
    char(di1(5)), char(di1(6)), char(di1(7))};
answer = inputdlg(prompt,'Parameters',[1 1 1 1 1 1 1],defAns);
di1={answer{1}, answer{2},answer{3},answer{4},answer{5},...
    answer{6},answer{7}};
vel=str2double(answer{1});
origvel=vel;    % velocity in cm/s
refreshrate=str2double(answer{2});
vel=floor(vres*vel/vscreen);%%pixel columns/s
%%pixel columns/refresh interval ---make sure at least 1

if vel==0
    dist=0;
else
    dist=max(floor(vel/refreshrate), 1);%drifting speed
end
cycledur=str2double(answer{3});%duration of stimulus,default 30 sec
dotnum=str2double(answer{4})/2; %dot number for full screen
dotrad=str2double(answer{5}(1,:)); %radius of dots
XenStage=str2double(answer{6}); % Xenopus Laevis stage
bardur=str2double(answer{7}); % time of running bar stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct / load stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch sel
    
    case 1 %for dot stimulus
        hdotstim = hres;
        vdotstim = vres/2; % half screen
%         tstart=[];
        
    case 2 % for full screen moving bar stimulus
        barunit=vres/(vscreen*0.25);   %pixel columns
        barunit=max(floor(barunit/2)*2, 2); %round, even, at least 2
        allvres=ceil(vres*1/barunit)*barunit;
        blockunit=[ones(hres, barunit/2) zeros(hres, barunit/2)];
        blockstim=repmat(blockunit, [1, allvres/barunit])*255;
        blockstim=uint8(blockstim);
                
    case 3 % for single bar stimulus
        barunit=vres/(vscreen*0.25);   %pixel columns
        barunit=max(floor(barunit/2)*2, 2); %round, even, at least 2
        allvres=ceil(vres*1/barunit)*barunit;
        blockunit=[ones(hres, barunit/2) zeros(hres, barunit/2)];
        blockstim=repmat(blockunit, [1, allvres/barunit])*255;
        blockstim (1:600,1:600)=0;
        blockstim=uint8(blockstim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imview(stim)
%run stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
%  window=Screen('OpenWindow',0,0,[0,0,800,600]);
window=Screen('OpenWindow',2,0); % open screen    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('started\n') 
%dot stimulus begin 
n=length(dotrad);
% Start the acquisition.  
try
 switch sel
        case 1
            %dot stimulus begin 
            n=length(dotrad);
            % Start the acquisition.
            start(vidobj)  
            for j=1:n
                tic
                dotradius=dotrad(j);
                dotradius=floor(dotradius*vres/vscreen); %%vertical reslution(pixle)/vertical size (cm)
                cll=2*dotradius+3; %diameter of dots why adding 3?
                nmcll=floor(hdotstim/cll)*floor(vdotstim/cll);%number of dots on half screen
                if dotnum>nmcll
                    error ('number of dots is too large')
                end
                pos=randperm(nmcll);
                h2=floor(hdotstim/cll);
                v2=floor(vdotstim/cll);
                %for poscnt=1:dotnum
                %    elm=pos(poscnt);
                %    pos(pos==(elm+1)) = [];
                %    pos(pos==(elm-1)) = [];
                %end
                %No dots vertically adjacent
                hOffset = floor((hdotstim - h2*cll)/h2);
                vOffset = floor((vdotstim - v2*cll)/v2);
                fprintf('started\n') 
                matunit=zeros(h2,v2);
                for cnt=1:dotnum
                    elm=pos(cnt);
                    rw=ceil(elm/v2);
                    cl=elm-(rw-1)*v2;
                    matunit(rw,cl)=1;
                end
                bh=repmat((-dotradius-1:dotradius+1)', [1, cll]);
                bv=repmat((-dotradius-1:dotradius+1), [cll, 1]);
                dotunit=(bh.^2+bv.^2).^0.5<=dotradius; %unit dot stimulus
                nrow=size(dotunit, 1);
                %number of columns in the imagedmat
                ncol=size(dotunit, 2); %REVERT!
                %number of rows in the imagemat
                dotunit_wOffset=zeros(nrow + vOffset,ncol + hOffset);
                dotunit_wOffset(1:nrow,1:ncol)=dotunit;
                dotunit=dotunit_wOffset;
                nrow=size(dotunit, 1);
                ncol=size(dotunit, 2);
                rSum = sum(matunit,2);
                for ir=1:size(rSum,1)
                    if rSum(ir)==0
                        matunit(ir,randi([1 size(matunit,2)]))=1;
                    end
                end
                cSum = sum(matunit,1);
                for ic=1:size(cSum,2)
                    if cSum(ic)==0
                        matunit(randi([1 size(matunit,1)]),ic)=1;
                    end
                end
                %check if any row/column has no dots
                stimhalf0=repmat(dotunit, size(matunit, 1), size(matunit, 2));
                [r c]=find(~matunit);
                for ii=1:length(r)
                    stimhalf0(nrow*(r(ii)-1) + 1 : nrow*(r(ii)-1) + nrow, ...
                        ncol*(c(ii)-1) + 1 : ncol*(c(ii)-1) + ncol)=0;
                end                
                if size(stimhalf0, 1)>600
                    stimhalf0=stimhalf0(1:600,1:end);
                end                
                hres=size(stimhalf0,1);
                vres=size(stimhalf0,2);
                stimhalf1=zeros(hres,vres);
                stimhalf1(1:hres,1:vres)=stimhalf0;
                stimhalf1=stimhalf1*255;
                stimfull=zeros(hres, vres*2);  %matrix for full screen display       
                stimfull(1:hres, 1:vres,:)=stimhalf1; % matrix for left drifting
                stimfull(1:hres, 401+400-vres:800,:)=stimhalf1;% matrix for right drifting
                stimhalf1=uint8(stimhalf1);
                stimfull=uint8(stimfull);
                sortid=num2str(j);%number of saved images
                pause(waitsec); %pauses 30 seconds before showing stimulus;
                t=cputime;                
                pos=1;
                c=0;
                toc
                tstart(1)=toc;  
                tic
                for i=t:t+cycledur*30 %duration of stimulus 30 sec
                    if pos+hres-1>hres
                        driftpos1=hres;
                        driftpos2=(pos+hres-1)-hres;
                    else driftpos1=pos+hres-1;
                        driftpos2=0;
                    end
                    if pos+dist>hres
                        nextpos=pos+dist-hres;
                    else
                        nextpos=pos+dist;
                    end
                    % half drifting image
                    stimdrift=[stimhalf1(pos:driftpos1,1:vres,:); ...
                        stimhalf1(1:driftpos2,1:vres,:)];
                    stimfull(1:hres,  401+400-vres:800)=stimdrift;
                    stimfull(1:hres,  1:vres)=stimdrift;
                    Screen(window,'PutImage',stimfull);
                    Screen('Flip', window);
                    pos=nextpos;
                    c=c+1;
                end
                toc
                tstart(2*j)=tstart(2*j-1)+toc;   
                %saving images before showing stimulus
                %end of saving images before showing stimulus
%                 Screen(window,'FillRect',0); %back to black background
%                 Screen('Flip', window);
%                 pause(switchinterval); % wait for 30s
                % switch drifting side
                tic
                for i=t:t+cycledur*30 %duration of stimulus 30 sec
                    if pos+hres-1>hres
                        driftpos1=hres;
                        driftpos2=(pos+hres-1)-hres;
                    else driftpos1=pos+hres-1;
                         driftpos2=0;
                    end
                    if pos+dist>hres
                        nextpos=pos+dist-hres;
                    else
                        nextpos=pos+dist;
                    end
                    % half drifting image
                    stimdrift=[stimhalf1(pos:driftpos1,1:vres,:); ...
                        stimhalf1(1:driftpos2,1:vres,:)];
                    stimfull(1:hres,  401+400-vres:800)=stimdrift;
                    stimfull(1:hres, 1:vres)=stimdrift;
                    Screen(window,'PutImage',stimfull);
                    Screen('Flip', window);
                    pos=nextpos;
                    c=c+1;
                end
                toc      
                tstart(2*j+1)=tstart(2*j)+toc; 
                %end of saving images before showing stimulus
                Screen(window,'FillRect',0); %back to black background
                Screen('Flip', window);
            end
            tstart(2*j+2)=tstart(2*j+1)+toc; 
            Screen('CloseAll');
            % save parameters to excel file            
            fid=fopen((['AvTest','_', char(dotorgrid(sel)),...
                 '.csv']),'wt');            
            x={'Dot Radius(cm)','tstart'};
            y={dotrad};
            z={tstart};
            [row1,col1]=size(x);[row2,col2]=size(y);[row3,col3]=size(z);
            for i=1:row1
                fprintf(fid,'%s,',x{i,1:end-1});
                fprintf(fid,'%s\n',x{i,end});
            end
            for i=1:row2
                fprintf(fid,'%6.4f,',y{i,1:end});
            end
            fprintf(fid,'\n');
            for i=1:row3
                fprintf(fid,'%6.4f,',z{i,1:end});
            end
            fclose(fid);
            % end in saving parameters to excel file
            %dot sitmulus end
        case 2
            % bar stimulus begin 
            start(vidobj)
            tic
            pos=1;
            c=0; 
            t=cputime;
            szv=size(blockstim,2);
            allszv=size(blockstim,2);
            pause(switchinterval); % wait for 30s
            toc
            tstart(1)=toc;
            tic
            for i=t:t+bardur*60 %duration of stimulus 30 sec
                if pos+hres-1>hres
                    driftpos1=hres;
                    driftpos2=(pos+hres-1)-hres;
                else driftpos1=pos+hres-1;
                    driftpos2=0;
                end
                if pos+dist>hres
                    nextpos=pos+dist-hres;
                else
                    nextpos=pos+dist;
                end
                % half drifting image
                blockstim1=[blockstim(1:end,pos:driftpos1,:) ...
                    blockstim(1:end,1:driftpos2,:)];
                Screen(window,'PutImage',blockstim1);
                Screen('Flip', window);
                pos=nextpos;
                c=c+1;                
            end                       
            Screen('CloseAll');
            toc 
            tstart(2)=tstart(1)+toc;
            % bar stimulus end
            % save parameters to excel file            
            fid=fopen((['AvTest','_', char(dotorgrid(sel)),...
                 '.csv']),'wt');            
            x={'tstart'};  z={tstart};
            [row1,col1]=size(x);[row3,col3]=size(z);
            for i=1:row1
                fprintf(fid,'%s,',x{i,1:end-1});
                fprintf(fid,'%s\n',x{i,end});
            end
            
            for i=1:row3
                fprintf(fid,'%6.4f,',z{i,1:end});
            end
            fclose(fid);
            % end in saving parameters to excel file
        case 3
            start(vidobj) 
            % single bar stimulus begin
            t=cputime;
            hres=size(blockstim,2);
            pos=1;
            c=0; 
            tic
            for i=t:t+bardur*30 %duration of stimulus 30 sec
                if pos+hres-1>hres
                    driftpos1=hres;
                    driftpos2=(pos+hres-1)-hres;
                else driftpos1=pos+hres-1;
                    driftpos2=0;
                end
                if pos+dist>hres
                    nextpos=pos+dist-hres;
                else
                    nextpos=pos+dist;
                end
                % half drifting image
                blockstim1=[blockstim(1:end,pos:driftpos1,:) ...
                    blockstim(1:end,1:driftpos2,:)];
                Screen(window,'PutImage',blockstim1);
                Screen('Flip', window);
                pos=nextpos;
                c=c+1;
            end
            toc         
            Screen('CloseAll');
            % bar stimulus end
    end
catch
    Screen('CloseAll');
    rethrow(lasterror);
    
end



% Wait for the acquisition to finish.
% wait(vidobj, 100)
% Once all frames have been written, close the file.
while (vidobj.FramesAcquired ~= vidobj.FramesPerTrigger) %vidobj.DiskLoggerFrameCount
    vidobj.FramesAcquired
    pause(5)
end

vidobj.FramesAcquired
vidobj.DiskLoggerFrameCount

close(vidobj.DiskLogger);
delete(vidobj)
clear vidobj
