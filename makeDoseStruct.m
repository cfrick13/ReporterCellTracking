function dosestruct = makeDoseStruct
close all
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
    exportdir = strcat(parentdir,'Tracking/Export');
cd(exportdir);

FileName = uigetfile('*export.mat');%choose file to load
% load(FileName)

%load metadata
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*Data.mat');
cd(exportdir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    load(char(filelist.name));
else
    filename = uigetfile();
    load(filename)
end


dosestruct = struct();
numberOfScenes = datastruct.sceneCount;
for i=1:numberOfScenes
   sceneStr = 's00';
   scene = num2str(i);
   sceneStr(end-(length(scene)-1):end)=scene;
   dosestruct(i).scene = sceneStr;
end


%input the number of doses
prompt = {'number of doses'};
dlg_title = 'enter the number of doses used in this experiment...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
numberOfDoses = str2num(inputdlgOutput{1});


%input the Tgfbeta concentrations of each dose
dosepromptarray = cell(1,numberOfDoses);
for i=1:numberOfDoses
    dosepromptarray{i} = strcat('dose',num2str(i));
end
prompt = dosepromptarray;
dlg_title = 'what are the doses?...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
doses = cellfun(@str2num,inputdlgOutput,'UniformOutput',1);







%input numbers for the scene numbers corresponding to each dose
doseandscenepromptarray = cell(1,numberOfDoses);
for i=1:numberOfDoses
    doseandscenepromptarray{i} = strcat('dose',num2str(i),'-',num2str(doses(i)));
end
prompt = doseandscenepromptarray;
dlg_title = 'what scenes correspond to which dose? (Use Matlab input notation "3,9:11")...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
doseToScene = cellfun(@str2num,inputdlgOutput,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the
%n doses

%input dose information into the structure
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
coloringArray = vertcat({dosestruct.(coloringChoice)});

for j=1:length(doseToScene)
    doseToScenemat = doseToScene{j};
    doseToSceneArray=[];
    for i = 1:length(doseToScenemat)
        dosestr = num2str(doseToScenemat(i)); 
        if length(dosestr)>1
            doseToSceneArray{i} = strcat('s',dosestr);
        else
            doseToSceneArray{i} = strcat('s0',dosestr); 
        end
    end
    
    
    disp(doseToSceneArray)

indicesChoice =channelregexpmaker(doseToSceneArray);
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
disp(dmat)
indices = ~dmat;
[dosestruct(indices).dose] = deal(doses(j));
end





%%%%%%%%%%%%%%%%%%%%%%  number of Tgfbeta additions %%%%%%%%%%%%%%%%%%
prompt = 'how many times was tgf added?...';
dlg_title = 'how many times was tgf added?...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
tgftime = cellfun(@str2num,inputdlgOutput,'UniformOutput',1);

%input numbers for the frame(s) after which tgfbeta was added
tgftimeandscenearray = cell(1,tgftime);
for i=1:tgftime
    tgftimeandscenearray{i} = strcat('addition#',num2str(tgftime(i)));
end
prompt = tgftimeandscenearray;
dlg_title = 'after what frame(s) was tgf added...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
tgfFrame = cellfun(@str2num,inputdlgOutput,'UniformOutput',1);

%input numbers for the scenes for each tgfbeta addition
tgftimeandscenearray = cell(1,tgftime);
for i=1:tgftime
    tgftimeandscenearray{i} = strcat('addition#',num2str(tgftime(i)),'-',num2str(tgfFrame(i)));
end
prompt = tgftimeandscenearray;
dlg_title = 'what scenes correspond to each addition each addition?(Use Matlab input notation "3,9:11")...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
tgfScenes = cellfun(@str2num,inputdlgOutput,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the
%n doses

%input dose information into the structure
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
coloringArray = vertcat({dosestruct.(coloringChoice)});

for j=1:length(tgfScenes)
    tgfToScenemat = tgfScenes{j};
    tgfToSceneArray=[];
    for i = 1:length(tgfToScenemat)
        dosestr = num2str(tgfToScenemat(i)); 
        if length(dosestr)>1
            tgfToSceneArray{i} = strcat('s',dosestr);
        else
            tgfToSceneArray{i} = strcat('s0',dosestr); 
        end
    end
    
    

indicesChoice =channelregexpmaker(tgfToSceneArray);
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
disp(dmat)
indices = ~dmat;
[dosestruct(indices).tgfFrame] = deal(tgfFrame(j));
end



%%%%%%%%%%%%%%%%%%%%%%  reference frames %%%%%%%%%%%%%%%%%%
prompt = 'what are the flatfield reference scenes?...';
dlg_title = 'what are the flatfield reference scenes?...';
inputdlgOutput = inputdlg(prompt,dlg_title,[1,80]);
BACKGROUND = cellfun(@str2num,inputdlgOutput,'UniformOutput',0);


for j=1:length(BACKGROUND)
    tgfToScenemat = BACKGROUND{j};
    tgfToSceneArray=[];
    for i = 1:length(tgfToScenemat)
        dosestr = num2str(tgfToScenemat(i)); 
        if length(dosestr)>1
            tgfToSceneArray{i} = strcat('s',dosestr);
        else
            tgfToSceneArray{i} = strcat('s0',dosestr); 
        end
    end
    
    

indicesChoice =channelregexpmaker(tgfToSceneArray);
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
disp(dmat)
indices = ~dmat;
[dosestruct(indices).flatfield] = deal('BACKGROUND');
[dosestruct(dmat).flatfield] = deal('experiment');
end









[~,b] = regexp(FileName,'exp[0-9]');
savename = FileName(1:b);
cd(exportdir)
savename = strcat(savename,'-DoseAndSceneData.mat');
save(savename)


    
    
    
    
    
    



stophere=1;


end





function channelinputs =channelregexpmaker(channelstoinput)
    channelinputs = '(';
    for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
        if i ==1
        channelinputs = strcat(channelinputs,channelstoinput{i});
        elseif i < length(channelstoinput)
            channelinputs = strcat(channelinputs,'|',channelstoinput{i});
        else
            channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
        end
    end
    
    if strcmp(')',channelinputs(end))
    else
        channelinputs(end+1)=')';
    end
end




