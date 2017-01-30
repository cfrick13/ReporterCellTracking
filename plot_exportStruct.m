function plot_exportStruct
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
cd(strcat(parentdir,'Tracking\Export'));

A = uigetfile;%choose file to load
load(A)

timeInterval = 4; %minutes
stimulationFrame = 1; %frame immediately preceeding stimulation
cellTracesString = 'medianNucEGFP'; %value to plot
finalFrame = 50; %choose which you want to be your final frame


    
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
% indicesChoice = {'s06','s11','s19'};
% indicesChoice = 's06|s11|s19';
indicesChoice = 's010';
coloringArray = vertcat({exportStruct.(coloringChoice)});
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
    uniqueColoring = unique(coloringArrayTrunc);
    cmap = lines(length(uniqueColoring));
    cmap(1,:) = [0 0 0];
    cmap(2,:) = [0 0.5 0.2];
    cmap(3,:) = [0.5 0.0 0.2];
disp(cmap)

colorNames = {'k','r','b','g'};
colormapMatrix = zeros(length(coloringArrayTrunc),size(cmap,2));
for i=1:length(coloringArrayTrunc)
   cA = coloringArrayTrunc{i};
   idx = strcmp(uniqueColoring,cA);
   colormapMatrix(i,:) = cmap(idx,:);
   colormapArray{i} = colorNames{idx};
   
end


%extract the cell traces for the desired number of frames
cellTracesFull = vertcat(exportStruct(indices).(cellTracesString));
cellTraces = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]

%normalize by basal values
basalVector = cellTraces(:,stimulationFrame); %88x1;
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
cellTracesNorm = cellTraces.*(invBasalMatrix');

%build the time vector and matrix for plotting
timeVector = (linspace(0,size(cellTraces,2)-1,size(cellTraces,2))-stimulationFrame).*timeInterval;
timeOnes = ones(size(cellTraces,1),1);
timeMatrix = timeOnes*timeVector;


f = figure(3)
    subplot(2,2,1);
        p = plot(timeMatrix',cellTraces','LineWidth',1);
            set(p, {'color'}, num2cell(colormapMatrix,2));
            xlim([timeVector(1) timeVector(end)])
            xlabel('minutes');
            ylabel(cellTracesString)
    subplot(2,2,2);
        p = plot(timeMatrix',cellTracesNorm','LineWidth',1);
            set(p, {'color'}, num2cell(colormapMatrix,2));
            xlim([timeVector(1) timeVector(end)])
            ylim([0 10])
            xlabel('minutes');
            ylabel(strcat('fold change'))
f.Position = [182 129 1058 969];
        




end