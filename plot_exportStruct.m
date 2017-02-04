function plot_exportStruct
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
load(FileName)

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


timeInterval = 4; %minutes
stimulationFrame = 10; %frame immediately preceeding stimulation
smadTracesString = 'medianNucEGFP'; %value to plot
reporterTracesString = 'totalNucRFP';

numberOfFrames = size(timeVec,2);
finalFrame = size(timeVec,2); 
finalFrame = 40; %choose which you want to be your final frame


    
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
indicesChoiceArray = {'s01','s02','s03','s04','s05','s06'};
indicesChoice =channelregexpmaker(indicesChoiceArray);

% indicesChoice = 's06|s11|s19';
% indicesChoice = 's010';
coloringArray = vertcat({exportStruct.(coloringChoice)});
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
    uniqueColoring = unique(coloringArrayTrunc);
%     cmap = colormap(viridis(length(uniqueColoring)));
    cmap = colormap(parula(length(uniqueColoring).*2)./2);
%     cmap(1,:) = [0 0 0];
%     cmap(2,:) = [0 0.5 0.2];
%     cmap(3,:) = [0.5 0.0 0.2];
disp(cmap)



colorNames = {'k','r','b','g','c'};
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
setTequalZeroToStimulation = timeVector(stimulationFrame);
xtickTimeVector = timeVector - setTequalZeroToStimulation;
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);
for i=1:length(coloringArrayTrunc)
   cA = coloringArrayTrunc{i};
   idx = strcmp(uniqueColoring,cA);
   timeMatrix(i,:) = timeVec(idx,1:finalFrame)-setTequalZeroToStimulation;  
end



%function to exract the cell traces, normalized and not
[smadCellTracesNorm,smadCellTraces] = extractTraces(exportStruct,indices,smadTracesString,finalFrame,stimulationFrame);
[reporterCellTracesNorm,reporterCellTraces] = extractTraces(exportStruct,indices,reporterTracesString,finalFrame,stimulationFrame);

%build the time vector and matrix for plotting
% timeVector = (linspace(0,size(smadCellTraces,2)-1,size(smadCellTraces,2))-stimulationFrame).*timeInterval;
% timeVector =1;
% timeOnes = ones(size(smadCellTraces,1),1);
% timeMatrix = timeOnes*timeVector;


f = figure(3);
    subplot(2,2,1);
        p = plot(timeMatrix',smadCellTraces','LineWidth',1.5,'Color',[0. 0 0]);
%             set(p, {'color'}, num2cell(colormapMatrix,2));
            xlim([-30 200])
            xlabel('minutes');
            ylabel('total nuclear fluorescence (au)')
            title('Level of endogenous nuclear NG-Smad3');
            h=p.Parent;
            h.XTick = [-30:30:400];
    subplot(2,2,2);
        p = plot(timeMatrix',smadCellTracesNorm','LineWidth',1.5,'Color',[0 0. 0]);
%             set(p, {'color'}, num2cell(colormapMatrix,2));
            xlim([-30 200])
            ylim([0 5])
            xlabel('minutes');
            ylabel(strcat('fold-change in total nuclear fluorescence'))
            h=p.Parent;
            h.XTick = [-30:30:400];
            title('Fold-change of endogenous nuclear NG-Smad3')
            t = text(0,0,strcat('N=',num2str(numberOfCells)));
            t.Units='normalized';
            t.Position = [0.95 0.9];
            t.HorizontalAlignment = 'right';
            t.FontSize = 12;
                t = text(0,0,strcat('2.5 ng/mL Tgf'));
                t.Units='normalized';
                t.Position = [0.95 0.85];
                t.HorizontalAlignment = 'right';
                t.FontSize = 12;
            
            
%     subplot(2,2,3);
%     scatter(smadCellTraces(:,1),reporterCellTraces(:,1));
%     subplot(2,2,4);
%     scatter(smadCellTraces(:,40),reporterCellTraces(:,40));
f.Position = [182 129 1400 1100];

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
basalVector = cellTraces(:,stimulationFrame); %88x1;
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
cellTracesNorm = cellTraces.*(invBasalMatrix');
end