% This is for batch running multiple movies with the TAD9000 function
clear

directory = uigetdir;
cd(directory);

mov_list = dir('*avi');

for i = 1:length(mov_list)
    [~,mov_name,~] = fileparts(mov_list(i).name);
    mov = VideoReader(mov_list(i).name);
    
    mkdir(mov_name)
    cd(mov_name)
    
    try
        [encAvg{i},numEncount{i},numAvoid{i}] = TadFunctionTest(mov);
        
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