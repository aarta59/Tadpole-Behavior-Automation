 % This is for batch running multiple movies with the TAD9000 function
clear

directory = uigetdir;
cd(directory);

mov_list = dir('*avi');

def = {'60','4.5','0.88','20','10','10','15','70','110'};
prompt = {'Gaussian hsize','Gaussian standard deviation (sigma)',...
    'Detection thresholding value', 'Tadpole speed variability',...
    'Tadpole X direction noise', 'Tadpole Y direction noise',...
    'Tadpole distance between eyes and gut', 'Lower avoidance angle tolerance',...
    'Higher avoidance angle tolerance'};
defAns = {char(def(1)), char(def(2)), char(def(3)), char(def(4)),...
    char(def(5)), char(def(6)), char(def(7)), char(def(8)), char(def(9))};
answer = inputdlg(prompt,'Initial Values', [1 1 1 1 1 1 1 1 1], defAns);
answer = str2double(answer);

for i = 1:length(mov_list)
    [~,mov_name,~] = fileparts(mov_list(i).name);
    mov = VideoReader(mov_list(i).name);
    
    mkdir(mov_name)
    cd(mov_name)
    
    try
        [encAvg{i},numEncount{i},numAvoid{i}] = TadFunctionTest(mov,answer);
        
        f_name = string(zeros(1,length(encAvg{i})));
        for j = 1:length(encAvg{i})
            f_name(1,j) = string(mov_name);
        end
    
        v_name{i} = f_name;
    
        cd(mov.path)
    catch e
        fprintf(2,'An error occurred in processing %s.\n',mov_name)
        fprintf(2,'The error message was:\n%s\n',e.message)
        cd(mov.path)
        continue
    end
    
end

name = cat(2,v_name{1:end});
avg = cat(2,encAvg{1:end});
enc = cat(2,numEncount{1:end});
avo = cat(2,numAvoid{1:end});

alldata = [name',avg',avo',enc'];
    
    
tab = table(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,4),'VariableNames',...
     {'MovieName','AvoidanceIndex','NumberAvoidances', 'NumberEncounters'});

writetable(tab,'Batch_Avoidance_Data','FileType','spreadsheet')
