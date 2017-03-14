function SigFeaturesArray = convert_exportStruct_doses_toInformation
close all
%determine the location of the matlab function and establish export
%directory in relation to that filepath
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
    exportdir = strcat(parentdir,'Tracking/Export');
cd(exportdir);

exportdirz = exportdir;
%load the exported tracking structure
FileName = uigetfile('*export.mat');%choose file to load
load(FileName)


%load metadata associated with the experiment (requires manual input if
%there is ambiguity
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*metaData.mat');
cd(exportdir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    filename = uigetfile();
    metaData = load(filename);
end

timeVec = metaData.timeVec;

%load information regarding doses and scenes and tgfbeta addition
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*DoseAndScene*');
cd(exportdirz)
filelist = dir(datequery);
    if isempty(filelist)
       dosestruct = makeDoseStruct; %run function to make doseStruct 
    else
        dosestructstruct = load(char(filelist.name));
        dosestruct = dosestructstruct.dosestruct;
    end
    
    %exportstruct
    %datastruct
    %dosestruct
    
    
    
    
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored 
colormapChoice = 'lines';
darkenFactor = 1.5;
    
%determine the scenes present in the experiment   
scenestr = 'scene';
sceneListArray = vertcat({exportStruct.(scenestr)});
sceneList = unique(sceneListArray);
sceneListArrayTwo = vertcat({dosestruct.(scenestr)});

%combine the exportStruct information with dosesstruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);


    dose = dosestruct(indicestwo).dose;
    frame = dosestruct(indicestwo).tgfFrame;
    
    dosestr = dosestruct(indicestwo).dosestr;
    framestr = dosestruct(indicestwo).tgfFramestr;
    
    
    [exportStruct(indices).dose] = deal(dose);
    [exportStruct(indices).frame] = deal(frame);
    [exportStruct(indices).dosestr] = deal(dosestr);
    [exportStruct(indices).framestr] = deal(framestr);
end

doseListArray = vertcat({exportStruct.dosestr});
doseList = unique(doseListArray);
    

    
%determine details needed for plotting such as when Tgfbeta is added, etc
%medianSmadbkg
stimulationFrame = exportStruct(1).frame;
smadTracesString = 'medianNucEGFP'; %value to plot
% smadTracesString = 'medianSmadbkg';
reporterTracesString = 'totalNucRFP';
numberOfFrames = size(timeVec,2);
finalFrame = numberOfFrames;




%establish the color map for plotting
coloringArray = vertcat({exportStruct.(coloringChoice)});
coloringList = unique(coloringArray);
indices = true(1,length(exportStruct));
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
    uniqueColoring = unique(coloringArrayTrunc);
    figure(1)
    cmap = colormap(parula(length(coloringList).*2));
    cmap = colormap(colormapChoice)./darkenFactor;
    close 1
    cmap = cmap; %darken the cmap


%assign a color array using the created colormap based on the choices above
colormapMatrix = zeros(length(coloringArrayTrunc),size(cmap,2));
for i=1:length(coloringArrayTrunc)
   cA = coloringArrayTrunc{i};
   idx = strcmp(uniqueColoring,cA);
   colormapMatrix(i,:) = cmap(idx,:);
%    colormapArray{i} = colorNames{idx};
   colormapArray{i} = cmap(idx,:);
end

%need to determine the number of scenes present and choose the time vector
%depending on the scene from which it was imaged
%THIS WORKS FOR NOW BUT NEEDS TO BE CHANGED
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);
    coloringArray = vertcat({exportStruct.(coloringChoice)});
    coloringList = unique(coloringArray);
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
for i=1:numberOfCells
    sceneChoice=exportStruct(i).scene;
    idx = strcmp(sceneListArray,sceneChoice);
    idxtwo = strcmp(sceneListArrayTwo,sceneChoice);
    
   stimulationFrame = dosestruct(idxtwo).tgfFrame;
   timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame);  
end
% setTequalZeroToStimulation = timeVector(stimulationFrame);
% xtickTimeVector = timeVector - setTequalZeroToStimulation;



indices = true(1,length(exportStruct));
%function to exract the cell traces, normalized and not
[smadCellTracesNorm,smadCellTraces] = extractTraces(exportStruct,indices,smadTracesString,finalFrame,stimulationFrame);
[reporterCellTracesNorm,reporterCellTraces] = extractTraces(exportStruct,indices,reporterTracesString,finalFrame,stimulationFrame);

%build the time vector and matrix for plotting
% timeVector = (linspace(0,size(smadCellTraces,2)-1,size(smadCellTraces,2))-stimulationFrame).*timeInterval;
% timeVector =1;
% timeOnes = ones(size(smadCellTraces,1),1);
% timeMatrix = timeOnes*timeVector;


response = struct();
SigFeaturesArray = cell(1,length(doseList));
for i=1:length(doseList)
    doseChoice = doseList{i};
    idx = strcmp(doseListArray,doseChoice);
    
    response.abslevel = smadCellTraces(idx,:);
    response.fc = smadCellTracesNorm(idx,:);
    
    fnames = fieldnames(response);
    responseMatrix = [];
    for jk=1:length(fnames)
        fn = fnames{jk};
        responseMatrix = horzcat(responseMatrix,response.(fn));
    end
        
     % SigFeaturesArray{1,i} = [n,m] where each m is a response type (measure) for n number of cells (observation)
%     responseMatrix = horzcat(abslevel,fc); %response matrix = [n,m] where each m is a response type (measure) for n number of cells (observation)
    SigFeaturesArray{1,i} = responseMatrix';% SigFeaturesArray = {1,n} where n = conditions (e.g. doses)
   
end
            

        




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
end


function [cellTracesNorm,cellTraces] = extractTraces(exportStruct,indices,xTracesString,finalFrame,stimulationFrame)
% extract the cell traces for the desired number of frames
cellTracesFull = vertcat(exportStruct(indices).(xTracesString));
cellTraces = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]

% normalize by basal values
basalVector = nanmedian(cellTraces(:,stimulationFrame-3:stimulationFrame),2);
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
cellTracesNorm = cellTraces.*(invBasalMatrix');
end