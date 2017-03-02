function plot_exportStruct_ReporterMedianTEST
close all
%determine the location of the matlab function and establish export
%directory in relation to that filepath
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
    loaddir = strcat(parentdir,'Export');
    exportdir = strcat(parentdir,'LookingAtData');
cd(loaddir);

exportdirz = exportdir;
%load the exported tracking structure
FileName = uigetfile('*export.mat');%choose file to load
load(FileName)



%load metadata associated with the experiment (requires manual input if
%there is ambiguity
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*metaData.mat');
figureSaveName = FileName(1:a-1);
cd(loaddir)
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
cd(loaddir)
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
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);
    coloringArray = vertcat({exportStruct.(coloringChoice)});
    coloringList = unique(coloringArray);
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
parfor i=1:numberOfCells
    sceneChoice=exportStruct(i).scene;
    idx = strcmp(sceneListArray,sceneChoice);
    idxtwo = strcmp(sceneListArrayTwo,sceneChoice);
    
   stimulationFramez = dosestruct(idxtwo).tgfFrame;
   timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFramez);  
end
% setTequalZeroToStimulation = timeVector(stimulationFrame);
% xtickTimeVector = timeVector - setTequalZeroToStimulation;



indices = true(1,length(exportStruct));
%function to exract the cell traces, normalized and not
[smadCellTracesNorm,smadCellTraces] = extractTraces(exportStruct,indices,smadTracesString,finalFrame,stimulationFrame);
[reporterCellTracesNorm,reporterCellTraces] = extractTraces(exportStruct,indices,reporterTracesString,finalFrame,stimulationFrame);

% % build the time vector and matrix for plotting
% % timeVector = (linspace(0,size(smadCellTraces,2)-1,size(smadCellTraces,2))-stimulationFrame).*timeInterval;
% % timeVector =1;
% % timeOnes = ones(size(smadCellTraces,1),1);
% % timeMatrix = timeOnes*timeVector;


basalTime = -15;
timeDuration = timeMatrix(1,end);
timeDuration = 90;

sceneChoices = {'(s03|s04)';'(s01|s02)'};
cd(exportdir);
scC = load('sceneChoices.mat');
% load('2017_02_08 plate exp1-test-test.mat')



plottingMat = reporterCellTraces;
plottingMatNorm = reporterCellTracesNorm;

sceneChoicesFull = scC.sceneChoices;
sccvec = 1:2:length(sceneChoicesFull);
for jmi = 1:length(sccvec)-1
    scc = sccvec(jmi);
    sccplus = sccvec(jmi)+1;
    sceneChoices = sceneChoicesFull(scc:sccplus);
    f = figure(jmi);
    for i=1:length(sceneChoices)
        sceneChoice = sceneChoices{i};
        idx = regexp(sceneListArray,sceneChoice);
        idx = ~cellfun(@isempty,idx,'UniformOutput',1);
        subplot(2,length(sceneChoices),i);
            p = plot(timeMatrix(idx,:)',plottingMat(idx,:)','LineWidth',1.5,'Color',[1 0.7 0.7]);hold on
            timeVector = timeMatrix(1,:);
            medianVector = nanmedian(plottingMat(idx,:),1);
            plot(timeVector,medianVector,'LineWidth',2,'Color',[0 0 0])
    %             set(p, {'color'}, num2cell(colormapMatrix(idx,:),2));
    %             xlim([-30 200])
    %             xlim([basalTime timeDuration])
                xlabel('minutes');
                ylabel('total nuclear fluorescence (au)')
                title(sceneChoice);
                h=p.Parent;
    %             h.XTick = [-30:30:400];
                ylim([0 max(medianVector).*2])
    %             ylim([0 60000])
        subplot(2,length(sceneChoices),i+length(sceneChoices));
            p = plot(timeMatrix(idx,:)',plottingMatNorm(idx,:)','LineWidth',1.5,'Color',[0.7 1 0.7]);hold on
            timeVector = timeMatrix(1,:);
            medianVector = nanmedian(plottingMatNorm(idx,:),1);
            plot(timeVector,medianVector,'LineWidth',2,'Color',[0 0 0])
    %             set(p, {'color'}, num2cell(colormapMatrix(idx,:),2));
    %             xlim([-30 200])
    %             xlim([basalTime timeDuration])
                ylim([0 5])
                xlabel('minutes');
                ylabel(strcat('fold-change'))
                h=p.Parent;
    %             h.XTick = [-30:30:400];
                title(strcat('Fold-change',sceneChoice));
                t = text(0,0,strcat('N=',num2str(sum(idx))));
                t.Units='normalized';
                t.Position = [0.95 0.9];
                t.HorizontalAlignment = 'right';
                t.FontSize = 12;
                    t = text(0,0,strcat(sceneChoice));
                    t.Units='normalized';
                    t.Position = [0.95 0.85];
                    t.HorizontalAlignment = 'right';
                    t.FontSize = 12;
    end

    %     subplot(2,2,3);
    %     scatter(smadCellTraces(:,1),reporterCellTraces(:,1));
    %     subplot(2,2,4);
    %     scatter(smadCellTraces(:,40),reporterCellTraces(:,40));
    f.Position = [182 200 1000 600];
    displayFile = FileName;
    [a] = regexp(displayFile,'_');
    displayFile(a) = '-';
    f=gcf;
    f.Units = 'normalized';
    txtsring = displayFile;
    t = uicontrol('Style','text','String',txtsring);
        t.Units = 'Normalized';
        t.Position = [0.01 0.93  0.15 0.05];
         t.FontSize = 14;
        t.HorizontalAlignment = 'left';


    numberOfCells = num2str(size(smadCellTraces,1));
    txtsring = strcat('N=',numberOfCells);
    t = uicontrol('Style','text','String',txtsring);
        t.Units = 'Normalized';
        t.Position = [0.01 0.91  0.1 0.05];
         t.FontSize = 14;
         t.HorizontalAlignment = 'left';


    exprimentDuration = num2str(round(timeMatrix(1,end)./60,2,'significant'));
    txtsring = strcat(exprimentDuration,'-hours');
    t = uicontrol('Style','text','String',txtsring);
        t.Units = 'Normalized';
        t.Position = [0.01 0.89  0.1 0.05];
        t.FontSize = 14;
        t.HorizontalAlignment = 'left';

    for h =f.Children'
        h.FontSize = 14;
        h.FontName = 'helvetica';
        stophere=1;
    end

    
    cd(exportdir)
    cd(strcat(figureSaveName,'-figures'))
    f=gcf;
sc = sceneChoices{1};
scc = sceneChoices{2};
scenesForFile = strcat(sc(2:4),sc(6:8),'-',scc(2:4),scc(6:8));
    filename = strcat(figureSaveName,'-',scenesForFile,'.fig');
    savefig(f,filename);
    filename = strcat(figureSaveName,'-',scenesForFile,'.png');
    saveas(f,filename)
    close(f);
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
basalLength = 3;
if (stimulationFrame-basalLength)<1
    basalLength=0;
end

basalVector = nanmedian(cellTraces(:,stimulationFrame-basalLength:stimulationFrame),2);
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
cellTracesNorm = cellTraces.*(invBasalMatrix');
end