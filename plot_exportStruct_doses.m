function plot_exportStruct_doses
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


%load the exported tracking structure
FileName = uigetfile('*export.mat');%choose file to load
load(FileName)


%load metadata associated with the experiment (requires manual input if
%there is ambiguity
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


%load information regarding doses and scenes and tgfbeta addition
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*DoseAndScene*');
cd(exportdir)
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
    
    
    
 %choose which field based upon which each cell trace will get colored    
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
stimulationFrame = exportStruct(1).frame;
smadTracesString = 'totalsNucEGFP'; %value to plot
reporterTracesString = 'totalNucRFP';
numberOfFrames = size(timeVec,2);
finalFrame = numberOfFrames;


    
coloringChoice = 'dosestr'; %choose which field based upon which each cell trace will get colored

%establish the color map for plotting
coloringArray = vertcat({exportStruct.(coloringChoice)});
coloringList = unique(coloringArray);
indices = true(1,length(exportStruct));
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
    uniqueColoring = unique(coloringArrayTrunc);
    figure(1)
    cmap = colormap(parula(length(coloringList).*2));
    close 1
    cmap = cmap./2; %darken the cmap
disp(cmap)


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
timeVector = ceil(timeVec(1,:));
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




f = figure(3);
for i=1:length(doseList)
    doseChoice = doseList{i};
    idx = strcmp(doseListArray,doseChoice);
    subplot(2,length(doseList),i);
        p = plot(timeMatrix(idx,:)',smadCellTraces(idx,:)','LineWidth',1.5,'Color',[0. 0 0]);
            set(p, {'color'}, num2cell(colormapMatrix(idx,:),2));
%             xlim([-30 200])
            xlim([-30 120])
            xlabel('minutes');
            ylabel('total nuclear fluorescence (au)')
            title('Level of endogenous nuclear NG-Smad3');
            h=p.Parent;
            h.XTick = [-30:30:400];
            ylim([0 max(max(smadCellTraces))])
%             ylim([0 5000])
    subplot(2,length(doseList),i+length(doseList));
        p = plot(timeMatrix(idx,:)',smadCellTracesNorm(idx,:)','LineWidth',1.5,'Color',[0 0. 0]);
            set(p, {'color'}, num2cell(colormapMatrix(idx,:),2));
%             xlim([-30 200])
            xlim([-30 120])
            ylim([0 5])
            xlabel('minutes');
            ylabel(strcat('fold-change in total nuclear fluorescence'))
            h=p.Parent;
            h.XTick = [-30:30:400];
            title('Fold-change of endogenous nuclear NG-Smad3')
            t = text(0,0,strcat('N=',num2str(sum(idx))));
            t.Units='normalized';
            t.Position = [0.95 0.9];
            t.HorizontalAlignment = 'right';
            t.FontSize = 12;
                t = text(0,0,strcat(doseChoice,' ng/mL Tgf'));
                t.Units='normalized';
                t.Position = [0.95 0.85];
                t.HorizontalAlignment = 'right';
                t.FontSize = 12;
end
            
%     subplot(2,2,3);
%     scatter(smadCellTraces(:,1),reporterCellTraces(:,1));
%     subplot(2,2,4);
%     scatter(smadCellTraces(:,40),reporterCellTraces(:,40));
f.Position = [182 510 2253 719];

for h =f.Children'
    h.FontSize = 14;
    h.FontName = 'helvetica';
    stophere=1;
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