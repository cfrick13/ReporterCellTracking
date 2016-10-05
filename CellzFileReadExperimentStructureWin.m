function CellzFileReadExperimentStructureWin
basallength=6;
fs=7;
tshift=0;
smoothing = 42;
DATAdir = 'D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled';
ADir = 'D:\Users\zeiss\Documents\';
% DATAdir = '\Users\frick\Documents\MATLAB\AllImagingDataCompiled';

cd (DATAdir); mkdir('saveSCALARS')
% dataName = 'AllCellTimeTracesNewest.mat';
% % dataNameTwo = 'ExperimentStructureNoTimeAdjustMat.mat';

%%%%%%%%
dataNameTwo = 'ExperimentStructure.mat';
EmmiesDataFile = 'ExperimentStructure.mat';

chosenspecies = '594_bhlhe40';
%%%%%%%%
ijij=1;

load(dataNameTwo);

COLOUR = determineCOLOUR;

MEEDS = calcMediansForDates(ExperimentStructure);
DataStruct = MedianCorrection(MEEDS,ExperimentStructure);

ScalarStruct = struct();
GeneScalarStruct = struct();
ScalarStructOnly = struct();
tracesstruct = struct();
Ttracesstruct =struct();
Mmrnastruct = struct();

runEmAll = fieldnames(DataStruct);
coouunt=1;
for SingleStruct = runEmAll';

   
    TimeTraceAndDeets = DataStruct.(char(SingleStruct));
    Dose = char(TimeTraceAndDeets.dose);
    
%     if strcmp(TimeTraceAndDeets.Date,'2014_09_30')
        if strcmp(TimeTraceAndDeets.positon,'s03')
            sopt=1;
        end      
%     disp(strcat(TimeTraceAndDeets.dose,TimeTraceAndDeets.positon));
%     coouunt=coouunt+coouunt+1;
%     disp(coouunt);
%     end
    
    peaklength = 52;
    Comms = TimeTraceAndDeets.cellComments;
    
    %%%
%     disp(Comms);
    %%%
    
    if 1<sum(strcmp(Comms,'focus'))
    elseif strcmp(Dose,'empty')
    elseif strcmp(Dose,'Empty')
    else

        % disp(TimeTraceAndDeets)
        % disp(TimeTraceAndDeets.cellTraces) 
        %%%%%%%%%%
        [SCALARS,NAMO,COLORS,PlotStruct] = CellzThat(peaklength,smoothing,TimeTraceAndDeets,basallength,tshift);
                              %CellzThat(peaklength,SMOOTHING,timeCoursesetc,%define the peak length if off to 20 minutes)
        %%%%%%%%%%%

        traces=PlotStruct.FC;
        abstraces = PlotStruct.CONCENTRATION;
        timeVector = TimeTraceAndDeets.timeVector;
        Frame = TimeTraceAndDeets.frame;
        FixTime = TimeTraceAndDeets.fixTime;
        Date = TimeTraceAndDeets.Date;
        Position = TimeTraceAndDeets.positon;Position(1) ='p';
        cellNumbers = TimeTraceAndDeets.cellNumbers;
        cellScalar = SCALARS{length(SCALARS)};

        %dose colors
        if strcmp(Dose,'0dot00');d=1;
        elseif strcmp(Dose,'0dot01');d=2;
        elseif strcmp(Dose,'0dot02');d=3;
        elseif strcmp(Dose,'0dot03');d=4;
        elseif strcmp(Dose,'0dot04');d=5;
        elseif strcmp(Dose,'0dot07');d=6;
        elseif strcmp(Dose,'0dot15');d=7;
        elseif strcmp(Dose,'0dot60');d=8;
        elseif strcmp(Dose,'2dot40');d=9;
        else d=9;end


        %%%%%%%%%%%%%%%%
        f = figure(1);
        sh = subplot(2,1,1); h = plot(timeVector(Frame-fs:end)-timeVector(Frame+1)+tshift,traces(:,Frame-fs:end));
        ylim([0 8])
        xlim([-16 60])
        shh = subplot(2,2,3); hh = plot(timeVector(Frame-fs:end)-timeVector(Frame+1)+tshift,traces(:,Frame-fs:end),'Color',COLOUR{d});hold on
        ylim([0 8])
        xlim([-16 60])
        shhh = subplot(2,2,4); hhh = plot(timeVector(Frame-fs:end)-timeVector(Frame+1)+tshift,traces(:,Frame-fs:end),'Color',COLOUR{d});hold on
        ylim([0 8])
        xlim([-16 60])
        %%%%%%%%%%%%%%%


        title(sh,strcat(char(Date),'...',char(Dose),'...',char(Position)));

        %%%%%%%%%%%%%%%
        ff = figure(2);
        subplot(2,2,1);
        hth = plot(timeVector(Frame-fs:end)-timeVector(Frame+1)+tshift,abstraces(:,Frame-fs:end),'Color',COLOUR{d}); hold off
        ylim([0 8000])
        xlim([-16 60])
        title(strcat(char(Date),'...',char(Dose),'...',char(Position)));

        subplot(2,2,3);
        hthh = plot(timeVector(Frame-fs:end)-timeVector(Frame+1)+tshift,traces(:,Frame-fs:end),'Color',COLOUR{d}); hold off
        ylim([0 8])
        xlim([-16 60])
        %%%%%%%%%%%%%%%%%%


        %%%%%%% LOAD EMMIES %%%%%%%%%%%%%%
        
        cd(strcat(ADir,'MATLAB\AllImagingDataCompiled'))
        filename = strcat('emmies*',Date,'*',chosenspecies,'*',Position);
        AAA = load(EmmiesDataFile, filename);
            EmFileName = fieldnames(AAA);
            EmFileChar = char(EmFileName);

        if isempty(EmFileName)
        emmies = nan(size(SCALARS{length(SCALARS)}'));
        else
        emmies = getfield(AAA,char(EmFileName));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        if strcmp(Position,'p03')
            stope=1;
        end
        [CommentStruct,leg] = CommentsIntoLegend(TimeTraceAndDeets,cellNumbers,cellScalar,emmies);


        %use CommentStruct logicals to determine which data to plot and include
        setLegend(sh,leg,FixTime);
        editPlotLineFormatting(h,emmies,CommentStruct)
        editPlotLineFormattingSecond(hh,emmies,CommentStruct)

        editPlotLineFormattingSecond(hth,emmies,CommentStruct)
        editPlotLineFormattingSecond(hthh,emmies,CommentStruct)

        editPlotLineFormattingSecond(hhh,emmies,CommentStruct)

        drawnow
        stophere=1;


        if ~isempty(EmFileName)
        %     currdir = pwd;
        %     cd('\Users\frick\Documents\')
        %     saveas(f,char(EmFileName),'jpeg');
        %     savefig(f,char(EmFileName));
        %     
        %     saveas(ff,strcat('mrna',char(EmFileName)),'jpeg');
        %     savefig(ff,strcat('mrna',char(EmFileName)));
        %     
        %     cd(currdir) 
        %     disp('smelly')
        % if strcmp(Position,'p24')
        % stophere=1;
        % end
        end


        if strcmp(Dose,'0dot00')
            stophere=1;
        end
        set(hhh,'Color',[0.7 0.7 0.7])
        % hth.delete;
        % hthh.delete;


        %%%%choose which cell subsets to keep\remove
          scalarMat = cell2mat(SCALARS');
        emrs = reshape(emmies,[length(emmies) 1]);
        emmies = emrs;
                if size(scalarMat,2) == size(emmies',2)
                scalarMat(end+1,:) = emmies';
                else
                    stophere=1;
                    emmies = nan(size(scalarMat(end,:)));
                    scalarMat(end+1,:) = emmies;
                end

        %%%COMMENT OUT TO INCLUDE ALL CELLS IN SCALAR DATA)
        scalarMatSelectLogical = scalarMatDataExcluder(CommentStruct,...
            {'OverDrivers','Dimmers','Remove','DoubleNuc','Saturated'});
        % scalarMatSelectLogical = scalarMatDataExcluder(CommentStruct,{'Remove','DoubleNuc','Saturated'});
        scalarMatSelect = scalarMat(:,scalarMatSelectLogical);


        if sum(isnan(scalarMatSelect(end,:))) == length(scalarMatSelect(end,:))
        else
            figure(55)
            subplot(3,1,2);
        % x = scatter(scalarMatSelect(5,:),scalarMatSelect(end,:),'MarkerFaceColor','none','MarkerEdgeColor',COLOUR{d});hold on
        x = scatter(scalarMatSelect(207,:),scalarMatSelect(end,:),'MarkerFaceColor','none','MarkerEdgeColor',COLOUR{d});hold on

        title(strcat(Position,Date,'fc'));
        ylim([0 500])
        xlim([0.5 8])
        stophere=1;
        x.DisplayName = strcat(Position,Date,Dose)

            figure(55)
            subplot(3,1,1);
        % x = scatter(scalarMatSelect(3,:),scalarMatSelect(end,:),'MarkerFaceColor','none','MarkerEdgeColor',COLOUR{d});hold on
        x = scatter(scalarMatSelect(177,:),scalarMatSelect(end,:),'MarkerFaceColor','none','MarkerEdgeColor',COLOUR{d});hold on
        title(strcat(Position,Date,'abs'));
        ylim([0 500])

        stophere=1;
        x.DisplayName = strcat(Position,Date,Dose);
        end


        DoseNamesPresent = fieldnames(tracesstruct);


        % ttraces = traces(scalarMatSelectLogical,:);
        ttraces = abstraces(scalarMatSelectLogical,:);

        if strcmp(Date,'2015_03_25')
            if strcmp(Position,'p27')
            stophere=1;
            end
        end

        %%%determine if structure exists for that gene name
        if sum(strcmp(DoseNamesPresent,strcat('i',Dose,'i')))
            traccprev = tracesstruct.(strcat('i',Dose,'i'));
            trac = ttraces(:,Frame-fs:end);
            tracc = zeros(size(trac,1),400);
            tracc(:,1:size(trac,2)) = trac;   
            tracconc = vertcat(traccprev,tracc);
            tracesstruct.(strcat('i',Dose,'i')) = tracconc;
        else
            trac = ttraces(:,Frame-fs:end);
            tracc = zeros(size(trac,1),400);
            tracc(:,1:size(trac,2)) = trac;
            tracconc = tracc;
            tracesstruct.(strcat('i',Dose,'i')) = tracconc;
        end

        %%%relative responses
        Ttraces = traces(scalarMatSelectLogical,:);
        %%%determine if structure exists for that gene name
        if sum(strcmp(DoseNamesPresent,strcat('i',Dose,'i')))
            Ttraccprev = Ttracesstruct.(strcat('i',Dose,'i'));
            Ttrac = Ttraces(:,Frame-fs:end);
            Ttracc = zeros(size(Ttrac,1),400);
            Ttracc(:,1:size(Ttrac,2)) = Ttrac;   
            Ttracconc = vertcat(Ttraccprev,Ttracc);
            Ttracesstruct.(strcat('i',Dose,'i')) = Ttracconc;
        else
            Ttrac = Ttraces(:,Frame-fs:end);
            Ttracc = zeros(size(Ttrac,1),400);
            Ttracc(:,1:size(Ttrac,2)) = Ttrac;
            Ttracconc = Ttracc;
            Ttracesstruct.(strcat('i',Dose,'i')) = Ttracconc;

        end

        %%%mRNA
        Mmrna = scalarMatSelect(end,:);
        %%%determine if structure exists for that gene name
        if sum(strcmp(DoseNamesPresent,strcat('i',Dose,'i')))
            Mmrnaprev = Mmrnastruct.(strcat('i',Dose,'i'));
            Mmrnac = Mmrna;
            Mmrnacconc = horzcat(Mmrnaprev,Mmrnac);
            Mmrnastruct.(strcat('i',Dose,'i')) = Mmrnacconc;
        else
            Mmrnastruct.(strcat('i',Dose,'i')) = Mmrna;
        end





        %     abstracesstruct.(Dose) = abstraces(:,Frame-8:end);

        %%%%
        % scalarMatSelect = scalarMat;      


        %%%%ConcatenateScalars
        [a,b] = regexp(EmFileChar,'(594_snail|594_smad7|594_pai1|647_snail|647_smad7|647_pai1|594_bhlhe40)');
        GeneName = EmFileChar(a:b);

        GeneNamesPresent = fieldnames(GeneScalarStruct);

        %%%determine if structure exists for that gene name
        if sum(strcmp(GeneNamesPresent,strcat('i',GeneName,'i')))
        ScalarStruct = GeneScalarStruct.(strcat('i',GeneName,'i'));

                DoseNamesPresent = fieldnames(ScalarStruct);
                if sum(strcmp(DoseNamesPresent,strcat('i',Dose,'i'))); 
                scalarold = ScalarStruct.(strcat('i',Dose,'i'));
                scalarMatConc =scalarold;
                    if ~(size(scalarMatSelect,1) ==  size(scalarMatConc,1))
                   stophere=1;
                    end
                ScalarStruct.(strcat('i',Dose,'i')) = horzcat(scalarMatSelect,scalarMatConc);
                else
                ScalarStruct.(strcat('i',Dose,'i')) = scalarMatSelect;
                end
                GeneScalarStruct.(strcat('i',GeneName,'i')) = ScalarStruct;

        else
                ScalarStruct.(strcat('i',Dose,'i')) = scalarMatSelect;
                GeneScalarStruct.(strcat('i',GeneName,'i')) = ScalarStruct;
        end





                DoseNamesPresent = fieldnames(ScalarStructOnly);
                if sum(strcmp(DoseNamesPresent,strcat('i',Dose,'i'))); 
                scalarold = ScalarStructOnly.(strcat('i',Dose,'i'));
                scalarMatConc =scalarold;
                ScalarStructOnly.(strcat('i',Dose,'i')) = horzcat(scalarMatSelect,scalarMatConc);
                else
                ScalarStructOnly.(strcat('i',Dose,'i')) = scalarMatSelect;
                end







        % filenamez = strcat('SCALARS.mat');
        % save(filenamez,'SCALARS','NAMO','COLORS');

    end
stophere=1;
end
stophere=1;

% cd('\Users\frick\Documents\MATLAB\AllImagingDataCompiled')
cd(strcat(ADir,'MATLAB\AllImagingDataCompiled'));
save('ScalarStructure.mat','ScalarStruct','NAMO','COLORS');
save(strcat('GeneScalarStructure',chosenspecies,EmmiesDataFile(1:end-4),'.mat'),'GeneScalarStruct','NAMO','COLORS');
save('ScalarStructureOnly.mat','ScalarStructOnly','NAMO','COLORS');
save('tracesstructabs.mat','tracesstruct');
save('tracesstructfc.mat','Ttracesstruct');
save('mrnastruct.mat','Mmrnastruct');
cd(strcat(ADir,'MATLAB'))

end



function scalarMatSelect = scalarMatDataExcluder(CommentStruct,excluderArray) %only keep non-overdrivers
for j = 1:length(excluderArray)
excluderLogical = CommentStruct.(excluderArray{j});
if j ==1
    excluderLogicalCombo = excluderLogical; 
else
    excluderLogicalCombo = vertcat(excluderLogical,excluderLogicalCombo);
end
end

scalarMatSum = sum(excluderLogicalCombo,1);
scalarMatSelect = ~logical(scalarMatSum);
end


function COLOUR = determineCOLOUR;
orange = [1 0.5 0];blue = [0 0.8 0.8];red = [0.9 0 0];yellow = [0.7 0.7 0];green = [0 0.9 0];purple = [0.9 0 0.9];cyan = [0.0 0.4 0.0];grey = [0.4 0.4 0.4];black = [0 0 0];
COLOUR = {black, red, orange, purple, yellow, blue, cyan, green, [0.1 0.2 0.1]}; 
end





function MEEDS = calcMediansForDates(ExperimentStructure);
StructDates = cellstr(fieldnames(ExperimentStructure));

MEEDS = struct();
for SDate = StructDates'
    meeds=[];
sDate = char(SDate);
ExpStructDatePositions = getfield(ExperimentStructure,sDate);
StructDatePositions = fieldnames(ExpStructDatePositions);
for SDPosition = StructDatePositions' 
sdPosition = char(SDPosition);
ExpStructDatePositionDetails = getfield(ExpStructDatePositions,sdPosition);

cellTraces = ExpStructDatePositionDetails.cellTraces;
timeVector = ExpStructDatePositionDetails.timeVector;
frame = ExpStructDatePositionDetails.frame;

sansbkg = cellTraces(3:end,:);
sansbkg(sansbkg==0)=NaN;
cellTraces(3:end,:) = sansbkg;
% cellTraces(cellTraces==0) = NaN;
backsubbed = backgroundsubtraction(cellTraces);
timecourses = vertcat(timeVector,backsubbed);
meedians = nanmean(backsubbed(2:2:end,frame-7:frame),2);
%%%%%
% meedians = ones(size(meedians));
%%%%%

if isempty(meeds)
    meeds = meedians;
else
    meeds = vertcat(meeds,meedians);
end

end
MEEDS.(sDate) = nanmedian(meeds);
end
stophere=1;
end

function DataStruct = MedianCorrection(MEEDS,ExperimentStructure)
%%%scale the median of absolute uniformly using a coefficient
mc=0.001;

DataStruct=struct();
StructDates = cellstr(fieldnames(ExperimentStructure));
for SDate = StructDates'
sDate = char(SDate);
ExpStructDatePositions = getfield(ExperimentStructure,sDate);
ExpStructDateMedian = getfield(MEEDS,sDate);
StructDatePositions = fieldnames(ExpStructDatePositions);
for SDPosition = StructDatePositions' 
sdPosition = char(SDPosition);
ExpStructDatePositionDetails = getfield(ExpStructDatePositions,sdPosition);

cellTraces = ExpStructDatePositionDetails.cellTraces;
timeVector = ExpStructDatePositionDetails.timeVector;
frame = ExpStructDatePositionDetails.frame;
ExpStructDatePositionDetails.Date = sDate(2:end-1);


sansbkg = cellTraces(3:end,:);
sansbkg(sansbkg==0)=NaN;
cellTraces(3:end,:) = sansbkg;
% cellTraces(cellTraces==0) = NaN;
backsubbed = backgroundsubtraction(cellTraces);
for ii = 2:2:size(backsubbed)
    backsubbed(ii,:) = backsubbed(ii,:)./(ExpStructDateMedian.*mc);
end

timecourse = vertcat(timeVector,backsubbed);

ExpStructDatePositionDetails.timecourse = timecourse;
DataStruct.(strcat(sDate,sdPosition)) = ExpStructDatePositionDetails;
end
end
end


function [CommentStruct,leg] = CommentsIntoLegend(TimeTraceAndDeets,cellNumbers,cellScalar,emmies)
%noisy
%delayed\delay



CommentStruct=struct();
CellComments = TimeTraceAndDeets.cellComments';
cellComments = CellComments(1:2:end);
OverDrivers = strcmp(cellComments,'overdriver');
Dimmers = strcmp(cellComments,'dimmer');
Remove = strcmp(cellComments,'remove');
DoubleNuc = strcmp(cellComments,'doublenuc');
    NospaceCFP = strcmp(cellComments,'no cfp');
    NoCFP = strcmp(cellComments,'nocfp');
NoCfp = logical(NospaceCFP+NoCFP);
Saturated = strcmp(cellComments,'saturated');

CommentStruct.OverDrivers = OverDrivers;
CommentStruct.Dimmers = Dimmers;
CommentStruct.Remove = Remove;
CommentStruct.DoubleNuc = DoubleNuc;
CommentStruct.NoCfp = NoCfp;
CommentStruct.Saturated = Saturated;

legOver = repmat({''},size(OverDrivers,1),size(OverDrivers,2));
legDim = legOver;
legRem = legOver;
legDub = legOver;
legNoC = legOver;
legSat = legOver;

legOver(OverDrivers) = {'OVERDRIVER'};
legOver(~OverDrivers) = {''};
legDim(Dimmers) = {'DIMMER'};
legDim(~Dimmers) = {''};
legRem(Remove) = {'REMOVED'};
legRem(~Remove) = {''};
legDub(DoubleNuc) = {'DOUBLENUC'};
legDub(~DoubleNuc) = {''};
legNoC(NoCfp) = {'NOCFP'};
legNoC(~NoCfp) = {''};
legSat(Saturated) = {'SATURATED'};
legSat(~Saturated) = {''};


a = length(cellNumbers);
b = length(emmies);
c = a==b;
if ~c
acn = reshape(cellNumbers,[length(cellNumbers) 1]);
acs = reshape(cellScalar,[length(cellScalar) 1]);  
acem = reshape(emmies,[length(emmies) 1]); 
prim =horzcat(acn,acs,nan(size(acn)));    
% prim = horzcat(cellNumbers,cellScalar',nan(size(cellNumbers)));
else  
acn = reshape(cellNumbers,[length(cellNumbers) 1]);
acs = reshape(cellScalar,[length(cellScalar) 1]);  
acem = reshape(emmies,[length(emmies) 1]); 
% prim =horzcat(cellNumbers,cellScalar',emmies);
prim =horzcat(acn,acs,acem);
end

primo = strcat((num2str(prim)),'...',cellstr(legOver'),' ',cellstr(legDim'),' ',cellstr(legRem'),' ',cellstr(legDub'),' ',cellstr(legNoC'),' ',cellstr(legSat'));
leg = char(primo);
end

function setLegend(sh,leg,FixTime);
lege =legend(sh,num2str(leg));
if strcmp(FixTime,'2hr')
lege = legend(sh,'2hr');
else
lege =legend(sh,num2str(leg));
end
set(lege,'FontSize',13)
set(lege,'Location','northeastoutside')
end

function backsubbed = backgroundsubtraction(DATA)
%Background subtracted data vector
backsubbed = zeros(size(DATA,1)-2,size(DATA,2)); %3 rows lost lose time row, and two background value rows

    %subtract CFP background
    for i=1:2:(size(DATA,1)-2);
    backsubbed(i,:) = DATA(i+2,:)-DATA(1,:);
    i+1;
    end
    %subtract mNG background
    for j=2:2:(size(DATA,1)-2);
    backsubbed(j,:) = DATA(j+2,:)-DATA(2,:);
    j+1;
    end
end


function editPlotLineFormatting(h,emmies,CommentStruct)
hline = findobj(h, 'type', 'line');
set(hline(~isnan(emmies)),'LineWidth',3)
set(hline(CommentStruct.OverDrivers),'LineStyle',':')
set(hline(CommentStruct.Dimmers),'LineStyle','--')
set(hline(CommentStruct.Remove),'LineStyle','none','Marker','.');
set(hline(CommentStruct.DoubleNuc),'Marker','d');
set(hline(CommentStruct.NoCfp),'LineStyle','-.');
end



function editPlotLineFormattingSecond(h,emmies,CommentStruct)
hline = findobj(h, 'type', 'line');
set(hline(1:end),'LineWidth',0.5);
set(hline(CommentStruct.OverDrivers),'LineStyle','none')
set(hline(CommentStruct.Dimmers),'LineStyle','none')
set(hline(CommentStruct.Remove),'LineStyle','none');
set(hline(CommentStruct.DoubleNuc),'LineStyle','none');
set(hline(CommentStruct.Saturated),'LineStyle','none');
% set(hline(CommentStruct.NoCfp),'LineStyle','none');
end



%% CellzThat
function [SCALARS,NAMO,COLORS,PlotStruct] = CellzThat(peakLEngth,smoothi,DataStruct,basallength,tshift)
PlotStruct = struct();
FRame = DataStruct.frame;
TimeVEctor = DataStruct.timeVector;
DATA = DataStruct.timecourse;

                    %%
                    %define colors
                    colourz.orange = [1 0.5 0]; colourz.blue = [0 0.8 0.8]; colourz.red = [0.9 0 0]; colourz.yellow = [0.7 0.7 0]; colourz.green = [0 0.9 0]; colourz.green = [0.9 0 0.9]; colourz.cyan = [0.0 0.4 0.0]; colourz.grey = [0.4 0.4 0.4]; colourz.black = [0 0 0];
                    orange = colourz.orange; blue = colourz.blue; red = colourz.red; yellow = colourz.yellow; green = colourz.green; green = colourz.green; cyan = colourz.cyan; grey = colourz.grey; black = colourz.black;

%%
TimeINterval = TimeVEctor(1,2)-TimeVEctor(1,1);
%Define the length of the timecourse analysis
peaklength = round(peakLEngth./TimeINterval)+FRame-1; %convert minutes to frames


%Raw data matrix
DATAS = DATA(2:1:size(DATA,1),:);

%compile timecourses to plot
[DATUMS,DATUMSNAME,timecourse] = timecourses(DATA,FRame,basallength);
mNG = timecourse.mNG;
CFP = timecourse.CFP;
YC = timecourse.YC;


%%
%determine which correction to make (CFP or moving average, etc) 
CONCENTRATION = chooseresponse(timecourse,smoothi);


%calculate basal mNG and basalCFP
basal = nanmean(CONCENTRATION(:,FRame-basallength:1:FRame),2);
basalcfp = nanmean(CFP(:,FRame-basallength:1:FRame),2);

SCALARS{1}= basal;
NAMO{1} = 'Basal Level';
COLORS{1} = orange;

%calculate basalcolor
basalcolor = calcbasalcolor(CONCENTRATION,colourz,basal);

%calcuate DIFFERENCE basal vetor for subtraction to generate absolute difference
for i=1:1:size(CONCENTRATION,1);
bs = (ones(size(CONCENTRATION)));
basalsubtract(i,:)= bs(i,:).*basal(i);
end
DIFFERENCE = CONCENTRATION-basalsubtract;

%Calculate FC
for i=1:1:size(mNG,1);
FC(i,:) = CONCENTRATION(i,:)./basal(i);
end


%% Calculate dynamical features of the timecourses
%calculate time to peak and peak FC. and all gradients

%generate smoothed functions of CONCENTRATION,DIFF and FC
for i= 1:1:size(FC,1)
yyFC(i,:) = smooth(TimeVEctor,FC(i,:),0.03,'moving');
yyABS(i,:) = smooth(TimeVEctor,CONCENTRATION(i,:),0.03,'moving');
yyDIFF(i,:) = smooth(TimeVEctor,DIFFERENCE(i,:),0.03,'moving');
end
[gradFC,rateFC,tresponseFC,fresponseFC] = computegradient(yyFC,TimeVEctor,peaklength,FRame,TimeINterval,tshift);
[gradABS,rateABS,tresponseABS,fresponseABS] = computegradient(yyABS,TimeVEctor,peaklength,FRame,TimeINterval,tshift);
[gradDIFF,rateDIFF,tresponseDIFF,fresponseDIFF] = computegradient(yyDIFF,TimeVEctor,peaklength,FRame,TimeINterval,tshift);

SCALARS{2} = tresponseFC;
NAMO{2}= 'Response Time';
COLORS{2} = blue;


%% calculate the absolute peak and absolute difference
[peakFC,tpeakFC,timetopeakFC] =  callpeak(FC,peaklength,FRame,TimeINterval);
[peakABS,tpeakABS,timetopeakABS] = callpeak(CONCENTRATION,peaklength,FRame,TimeINterval);
[peakDIFF,tpeakDIFF,timetopeakDIFF] = callpeak(DIFFERENCE,peaklength,FRame,TimeINterval);
peakPERCENT = peakFC - ones(size(peakFC));
 
    SCALARS{5} = peakFC;
    NAMO{5}='Max Level/Basal';
    COLORS{5} = green;
SCALARS{3} = peakABS;
NAMO{3} = 'Max Level';
COLORS{3} = red;
    SCALARS{4} = peakDIFF;
    NAMO{4} = 'Max Level - Basal';
    COLORS{4}=yellow;
SCALARS{6} = peakPERCENT;
NAMO{6} = 'Max Percent Change';
COLORS{6} = green;


%% calculate average after TGFbeta addition
for i=1:1:size(CONCENTRATION,1);
meanCONCENTRATION(i,:) = nanmean(CONCENTRATION(i,FRame:1:peaklength));
meanDIFFERENCE(i,:) = nanmean(DIFFERENCE(i,FRame:1:peaklength));
meanFC(i,:) = nanmean(FC(i,FRame:1:peaklength));
meanPERCENT(i,:) =nanmean(FC(i,FRame:1:peaklength)-1);
end
SCALARS{7} = meanCONCENTRATION;
NAMO{7}='Mean Level';
COLORS{7} = red;
    SCALARS{8}=meanDIFFERENCE;
    NAMO{8}='Mean Level - Basal';
    COLORS{8} = yellow;
SCALARS{9}=meanFC;
NAMO{9}='Mean Level/Basal';
COLORS{9} = green;
    SCALARS{10}=meanPERCENT;
    NAMO{10}='Mean Percent Change';
    COLORS{10} = green;
    

%% calculate average after reaching peak
for i=1:1:size(FC,1)
meanCONCENTRATIONafterpeak(i,:) = nanmean(CONCENTRATION(i,tpeakFC(i):peaklength));
meanDIFFERENCEafterpeak(i,:) = nanmean(DIFFERENCE(i,tpeakFC(i):peaklength));
meanFCafterpeak(i,:) = nanmean(FC(i,tpeakFC(i):peaklength));
meanPERCENTafterpeak(i,:) = nanmean(FC(i,tpeakFC(i):peaklength))-1;
end

SCALARS{11} = meanCONCENTRATIONafterpeak;
NAMO{11} = 'Mean Level after Maximum';
COLORS{11} = red;
    SCALARS{12} = meanDIFFERENCEafterpeak;
    NAMO{12} = 'Mean Level - Basal after Maximum';
    COLORS{12} = yellow;
SCALARS{13} = meanFCafterpeak;
NAMO{13} = 'Mean Level/Basal after Maximum';
COLORS{13} = green;
    SCALARS{14} = meanPERCENTafterpeak;
    NAMO{14} = 'Mean Percent Change after Maximum';
    COLORS{14} = green;



%% Find AUC after TGFaddition
for i=1:1:size(mNG,1);
aucCONCENTRATION(i) = trapz(CONCENTRATION(i,FRame:1:peaklength));
aucDIFFERENCE(i)  = trapz(DIFFERENCE(i,FRame:1:peaklength));
aucFC(i) = trapz(FC(i,FRame:1:peaklength));
aucPERCENT(i) = trapz(FC(i,FRame:1:peaklength)-1);
end

SCALARS{18} = aucCONCENTRATION;
NAMO{18} = 'Integrated Level';
COLORS{18} = red;
    SCALARS{19} = aucDIFFERENCE;
    NAMO{19} = 'Integrated Level-Basal';
    COLORS{19} = yellow;
SCALARS{20} = aucFC;
NAMO{20} = 'Integrated Level/Basal';
COLORS{20} = green;
    SCALARS{21} = aucPERCENT;
    NAMO{21} = 'Integrated Percent Change';
    COLORS{21} = green;
    
SCALARS{15} = rateABS;
NAMO{15}= 'Maximum d/dt(Level)';
COLORS{15} = [0 0 0];

% SCALARS{17} = rateDIFF./basal;
% NAMO{17}= 'WEBrate';
SCALARS{16} = rateDIFF;
NAMO{16}= 'Maximum d/dt(Level-Basal)';
COLORS{16} = [0 0 0];

SCALARS{17} = rateFC;
NAMO{17} = 'Maximum d/dt(Level/Basal)';
COLORS{17} = [0 0 0];


%%
t=0;

% t=-2;
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeFifty(CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'TimeFiftyUp',red,TimeINterval,FRame,peakLEngth);

% t=-2;
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = ResponsePostTime(CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'(Level)afterResponse',red,TimeINterval,FRame,peakLEngth);


% t=-2;
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = ResponsePostTime(DIFFERENCE,SCALARS,NAMO,COLORS,sn,t,'(Level-Basal)afterResponse',red,TimeINterval,FRame,peakLEngth);


% t=-2;
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = ResponsePostTime(FC,SCALARS,NAMO,COLORS,sn,t,'(Level/Basal)afterResponse',red,TimeINterval,FRame,peakLEngth);


%t=-2
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponses(CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'Level',red,TimeINterval,FRame,peakLEngth);

% t=-2;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponses(DIFFERENCE,SCALARS,NAMO,COLORS,sn,t,'(Level-Basal)',yellow,TimeINterval,FRame,peakLEngth);

% t=-2;
sn = length(NAMO)+1;
% [SCALARS,NAMO,COLORS] = timeResponses(FC,SCALARS,NAMO,COLORS,sn,t,'FC',green,x,y);
[SCALARS,NAMO,COLORS] = timeResponses(FC,SCALARS,NAMO,COLORS,sn,t,'(Level/Basal)',green,TimeINterval,FRame,peakLEngth);

% t=-1;
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponsesAveraged(CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'Level',red,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponsesAveraged(DIFFERENCE,SCALARS,NAMO,COLORS,sn,t,'(Level-Basal)',yellow,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
% [SCALARS,NAMO,COLORS] = timeResponses(FC,SCALARS,NAMO,COLORS,sn,t,'FC',green,x,y);
[SCALARS,NAMO,COLORS] = timeResponsesAveraged(FC,SCALARS,NAMO,COLORS,sn,t,'(Level/Basal)',green,TimeINterval,FRame,peakLEngth);

% t=-1;
sn=length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponsesRate(gradABS,SCALARS,NAMO,COLORS,sn,t,'d/dt(Level)',red,TimeINterval,2,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponsesRate(gradDIFF,SCALARS,NAMO,COLORS,sn,t,'d/dt(Level-Basal)',yellow,TimeINterval,2,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS] = timeResponsesRate(gradFC,SCALARS,NAMO,COLORS,sn,t,'d/dt(Level/Basal)',green,TimeINterval,2,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=timeResponsesIntegrated(CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'Integrated Level',red,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=timeResponsesIntegrated(DIFFERENCE,SCALARS,NAMO,COLORS,sn,t,'Integrated Level-Basal',yellow,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=timeResponsesIntegrated(FC,SCALARS,NAMO,COLORS,sn,t,'Integrated Level/Basal',green,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=timeResponsesEndAveraged(TimeVEctor,CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'Level',red,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=timeResponsesEndAveraged(TimeVEctor,DIFFERENCE,SCALARS,NAMO,COLORS,sn,t,'Level-Basal',yellow,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=timeResponsesEndAveraged(TimeVEctor,FC,SCALARS,NAMO,COLORS,sn,t,'Level/Basal',green,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=maxVarOvertime(TimeVEctor,CONCENTRATION,SCALARS,NAMO,COLORS,sn,t,'Max Level t',red,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=maxVarOvertime(TimeVEctor,DIFFERENCE,SCALARS,NAMO,COLORS,sn,t,'Max Level-Basal t',yellow,TimeINterval,FRame,peakLEngth);

% t=-1;
sn = length(NAMO)+1;
[SCALARS,NAMO,COLORS]=maxVarOvertime(TimeVEctor,FC,SCALARS,NAMO,COLORS,sn,t,'Max Level/Basal t',green,TimeINterval,FRame,peakLEngth);

%% save scalars namo and colors files
for jj = 1:length(SCALARS)
lookclose = SCALARS{jj};
if size(lookclose,1) == 1
    lookagain = lookclose;
else
    lookagain = lookclose';
end
SCALARS{jj} = lookagain;
end


PlotStruct.CONCENTRATION = CONCENTRATION;
PlotStruct.DIFFERENCE = DIFFERENCE;
PlotStruct.FC = FC;
end






%%
function [DATUMS,DATUMSNAME,timecourse] = timecourses(DATA,FRame,basallength)   
mNGandCFP = DATA(2:end,:);
%mNG vector
mNG = mNGandCFP(2:2:size(mNGandCFP,1),:);
    DATUMS{2} = mNG;
    DATUMSNAME{2} = 'mNG';
    timecourse.mNG = mNG;
%CFP vector
CFP = mNGandCFP(1:2:size(mNGandCFP,1),:);
    DATUMS{3} = CFP;
    DATUMSNAME{3} = 'CFP';
    timecourse.CFP = CFP;

%normalized mNG vector

mNGbasals = nanmean(mNG(:,FRame-basallength:1:FRame),2);
for i = 1:1:size(mNG,1);
mNG_N(i,:) = mNG(i,:)./mNGbasals(i);
i+1;
end
    DATUMS{5} =mNG_N;
    DATUMSNAME{5} = 'mNG_N';

%normalized CFP vector
CFPbasals = nanmean(CFP(:,FRame-basallength:1:FRame),2);
for i = 1:1:size(CFP,1);
CFP_N(i,:) = CFP(i,:)./CFPbasals(i);
i+1;
end
    DATUMS{6} = CFP_N;
    DATUMSNAME{6} = 'CFP_N';
    
    %Y/C vector
YC = mNG./CFP_N;
    DATUMS{4} = YC;
    DATUMSNAME{4} = 'YC';
    timecourse.YC = YC;

%normalized YC vector
for i = 1:1:size(YC,1);
YCbasals(i) = nanmean(YC(i,1:FRame));
YC_N(i,:) = YC(i,:)./YCbasals(i);
i+1;
end
    DATUMS{7} = YC_N;
    DATUMSNAME{7} = 'YC_N';

end  



function response = chooseresponse(timecourse,smoothi)
mNG = timecourse.mNG;
CFP = timecourse.CFP;
YC = timecourse.YC;

area = [0;1];
response =mNG(:,1:1:size(mNG,2));

if smoothi == 1; 
    response = mNG(:,1:1:size(mNG,2));%mNG
elseif smoothi == 2;
    response = CFP(:,1:1:size(CFP,2));%CFP
elseif smoothi == 12;
    response = YC(:,1:1:size(YC,2)); %YC = mNG/CFP
elseif smoothi == 3                      %area corrected
        area = aa;
            for i = 1:1:size(response,1)
            basalarea(i) = nanmean(area(i,1:1:y));
            area(i,:) = area(i,:)/basalarea(i);
            end
        response =mNG(:,1:1:size(mNG,2)).*area(:,1:1:size(mNG,2));
        CFP = CFP(:,1:1:size(CFP,2)).*area(:,1:1:size(CFP,2));
            for i = 1:1:size(response,1)
            basalarea(i) = nanmean(area(i,1:1:y));
            end
elseif smoothi==4                        %4 is moving average
    for i = 1:size(mNG,1)       
        for j = 1:size(mNG,2)-2
    MeanSet = ([mNG(i,j) mNG(i,j+1) mNG(i,j+2)]);
    Summa = sum(~isnan(MeanSet));
    response(i,j) = nansum(MeanSet)./Summa;
        end
    end
    
elseif smoothi==42                   %42 is moving average of YC

% % %         for j = 1:size(YC,2)-2
% % %     MeanSet = ([YC(i,j) YC(i,j+1) YC(i,j+2)]);
% % %     Summa = sum(~isnan(MeanSet));
% % %     response(i,j) = nansum(MeanSet)./Summa;
% % %         end
% % %         for j = 2:size(YC,2)-1
% % %     MeanSet = ([YC(i,j-1) YC(i,j) YC(i,j+1)]);
% % %     Summa = sum(~isnan(MeanSet));
% % %     response(i,j) = nansum(MeanSet)./Summa;
% % %         end
% % %         for j = 3:size(YC,2)
% % %     MeanSet = ([YC(i,j-2) YC(i,j-1) YC(i,j)]);
    for i = 1:size(YC,1)
% % % % % % % % % % % % 
%         for j = 3:size(YC,2)
%      MeanSet = ([YC(i,j-2) YC(i,j-1) YC(i,j)]);
%     Summa = sum(~isnan(MeanSet));
%     response(i,j) = nansum(MeanSet)./Summa;
%         end
%         
        for j = 3:size(YC,2)-1
     MeanSet = ([YC(i,j-1) YC(i,j) YC(i,j+1)]);
    Summa = sum(~isnan(MeanSet));
    response(i,j) = nansum(MeanSet)./Summa;
        end
% % % % % % % % % % % 
% %       for j = 2:size(YC,2)-1
% %     MeanSet = ([YC(i,j-1) YC(i,j+0) YC(i,j+1)]);
% %     if isnan(MeanSet) == [1 0 0]
% %         MeanSet = ([YC(i,j+0) YC(i,j+1) YC(i,j+2)]);
% %     elseif isnan(MeanSet) == [1 1 0]
% %         MeanSet = ([YC(i,j+1) YC(i,j+2) YC(i,j+3)]);
% %     elseif isnan(MeanSet) == [0 0 1]
% %         MeanSet = ([YC(i,j-2) YC(i,j-1) YC(i,j+0)]);
% %     elseif isnan(MeanSet) == [0 1 1]
% %         MeanSet = ([YC(i,j-3) YC(i,j-2) YC(i,j-1)]);
% % %     end
% %     Summa = sum(~isnan(MeanSet));
% %     response(i,j) = nansum(MeanSet)./Summa;
% %       end 
    end
end
end












%% OLD FUNCTIONS















function datesTogether = compileDatesTogether(timecourses)
dosesTogether = setintialemptyTogether(timecourses,'dose');
datesTogether = setintialTogether(timecourses,'date',dosesTogether);
for i = 1:length(timecourses)
    dose = getfield(timecourses,{i},'dose');
    date = getfield(timecourses,{i},'date');
    frame = str2num(getfield(timecourses,{i},'frame'));
    DATA = timecourses(i).TRACES;
    currentdatecomp = getfield(datesTogether,strcat('i',date,'i'));
    currentcomp = getfield(currentdatecomp,strcat('i',dose,'i'));
    M = compiledatadoses(DATA,dose,currentcomp,frame);
    currentdatecomp = setfield(currentdatecomp,strcat('i',dose,'i'),M);
    datesTogether = setfield(datesTogether,strcat('i',date,'i'),currentdatecomp);
end
end

function dosesTogether = compileDosesTogether(timecourses)
dosesTogether = setintialemptyTogether(timecourses,'dose');
for i = 1:length(timecourses)
    dose = getfield(timecourses,{i},'dose');
    frame = str2num(getfield(timecourses,{i},'frame'));
    DATA = timecourses(i).TRACES;
    currentcomp = getfield(dosesTogether,strcat('i',dose,'i'));
    M = compiledatadoses(DATA,dose,currentcomp,frame);
    dosesTogether = setfield(dosesTogether,strcat('i',dose,'i'),M);
end
end

function timecoursesMedCorr = medianCorrectTimecourses(timecourses,medians)
for i = 1:length(timecourses)
    DATE = getfield(timecourses,{i},'date');
    med = getfield(medians,char(strcat('i',DATE,'i')));
    traces = getfield(timecourses,{i},'TRACES');
    traces(4:end,:) = traces(4:end,:)./med;
    timecourses = setfield(timecourses,{i},'TRACES',traces);
end
timecoursesMedCorr = timecourses;
end

function timecoursesBackSub = backSubTimecourses(timecourses)
timecoursesBackSub = timecourses;
for i = 1:length(timecoursesBackSub)
    DATE = getfield(timecoursesBackSub,{i},'date');
    traces = getfield(timecoursesBackSub,{i},'TRACES');
    backsubbed = backgroundsubtraction(traces);
    TRACES = zeros(size(traces));
    TRACES(1,:) = traces(1,:);
    TRACES(4:end,:) = backsubbed;
    timecoursesBackSub = setfield(timecoursesBackSub,{i},'TRACES',TRACES);
end
end

function datesTogetherMedCorr = medianCorrectDates(datesTogetherMedCorr,DATES,DOSES,medians)
for DATE = DATES
    datescomp = getfield(datesTogetherMedCorr,char(strcat('i',DATE,'i')));
    med = getfield(medians,char(strcat('i',DATE,'i')));
    for DOSE = DOSES
        dosesofdatecomp = getfield(datescomp,char(strcat('i',DOSE,'i')));
        if isempty(dosesofdatecomp)
        else
        dosesofdatecomp_medcorr = dosesofdatecomp./med;
        datescomp = setfield(datescomp,char(strcat('i',DOSE,'i')),dosesofdatecomp_medcorr);
        end
    end
    datesTogetherMedCorr = setfield(datesTogetherMedCorr,char(strcat('i',DATE,'i')),datescomp);
end
end

function medians = getDateMedians(datesTogether,DATES,DOSES)
medians = struct();
for DATE = DATES
    datescomp = getfield(datesTogether,char(strcat('i',DATE,'i')));
    medcomp = [];
    for DOSE = DOSES
        dosesofdatecomp = getfield(datescomp,char(strcat('i',DOSE,'i')));
        if isempty(dosesofdatecomp)
        else
            if isempty(medcomp)
        medcomp = dosesofdatecomp(4:end,8);
            else
        medcomp = vertcat(medcomp,dosesofdatecomp(4:end,8));
            end
        end
    end
    med = median(medcomp);
    medians = setfield(medians,char(strcat('i',DATE,'i')),med);
end
end

function FIELDS = getAllField(timecourses,FIELDstring)
j=1;
FIELDS = [];
testfield=[];
for i = 1:length(timecourses)
field = getfield(timecourses,{i},FIELDstring);

if length(FIELDS)>0
for k = 1:length(FIELDS)
    testfield(k) = mean(FIELDS{k} == field);
end
else
    FIELDS{j} = field;
    j=j+1;
end

if max(testfield) < 1;
FIELDS{j} = field; 
j=j+1;
else
end

fieldprev = field;
end

end

function together = setintialemptyTogether(timecourses,FIELDstring)
together = struct();
M=[];
for i = 1:length(timecourses)
field = getfield(timecourses,{i},FIELDstring);
together = setfield(together,strcat('i',field,'i'),M);
end
end

function together = setintialTogether(timecourses,FIELDstring,dosestructure)
together = struct();
M = dosestructure;
for i = 1:length(timecourses)
field = getfield(timecourses,{i},FIELDstring);
together = setfield(together,strcat('i',field,'i'),M);
end
end

function traces = getNames(DATAdir,dataName)
cd (DATAdir)
load(dataName);
clear DATAdir dataName
traces = whos;
% load(dataName,'snail*');
% load(dataName,'snail*');
end

function timecourses = getTimes(dataName,snail)
cycleM=1;
timecourses = struct;
for cfile = {snail.name}
    filestring = char(cfile);
    
    [a,b] = regexp(filestring,'201.......');
    date = filestring(a:b);
    [c,d] = regexp(filestring,'(smad3g|smad3f)');
    celltype = filestring(c:d);
    [e,f] = regexp(filestring,'p[0-9]+');
    position = filestring(e:f);
    [g,h] = regexp(filestring,'[0-9]dot[0-9]+');
    dose = filestring(g:h);
    [k,l] = regexp(filestring,'frame[0-9]+');
    frame = filestring(l-1:l);
    
        A = load(dataName,filestring);
        M = getfield(A,filestring);
        
        timecourses = setfield(timecourses,{cycleM},'dose',dose);
        timecourses = setfield(timecourses,{cycleM},'position',position);
        timecourses = setfield(timecourses,{cycleM},'date',date);
        timecourses = setfield(timecourses,{cycleM},'celltype',celltype);
        timecourses = setfield(timecourses,{cycleM},'frame',frame);
        timecourses = setfield(timecourses,{cycleM},'TRACES',M);
        cycleM=cycleM+1;
end

end



function M = compiledata(DATA,DOSE,currentcomp)
    backsubbed = backgroundsubtraction(DATA);
    if isempty(currentcomp)
    TIME = DATA(1,:);
    COMPILED(1,:) = TIME;
    COMPILED(2:3,:) = zeros(2,size(TIME,2));
    COMP = COMPILED;
    else
    COMP = currentcomp;
    end
    COMP = vertcat(COMP,backsubbed);
    M = COMP;
end


function M = compiledatadoses(DATA,DOSE,currentcomp,frame)
    matrixlength = 400;
    backsubwithextra = NaN(size(DATA,1)-3,matrixlength);
    backsubbed = backgroundsubtraction(DATA);
    
    
    frameshift = (frame-7);
    endshift = 1+size(DATA,2)-frameshift;
    
    backsubwithextra(:,1:endshift) = backsubbed(:,frameshift:end);
    
    if isempty(currentcomp)
    COMPILED = NaN(3,matrixlength);
    TIME(1,:) = DATA(1,frameshift:end);
    COMPILED(1,1:endshift) = TIME;
    COMPILED(2:3,1:endshift) = zeros(2,size(TIME,2));
    COMP = COMPILED;
    else
    COMP = currentcomp;
    end
    COMP = vertcat(COMP,backsubwithextra);
    M = COMP;
end    

%%


%%
function h = plotAnalysisMatrices(TIME,PLOT,CONCENTRATION,basalcolor,FC,x,y,doseNumber)

frameshift = y-7;
TimeVector = TIME(frameshift:end) - ones(size(TIME(frameshift:end))).*min(TIME(frameshift:end));

% 
% figure(1)
% subplot(2,6,doseNumber);plot(TimeVector,CONCENTRATION(:,frameshift:end),'Color',basalcolor,'LineWidth',0.1);hold on
% ylim([0 8])
% xlim([0 100])
% subplot(2,6,6+doseNumber);plot(TimeVector,FC(:,frameshift:end),'Color',basalcolor,'LineWidth',0.1);hold on
% ylim([0 8])
% xlim([0 100])
% 
% 
% figure(2)
% subplot(2,1,1);plot(TimeVector,CONCENTRATION(:,frameshift:end),'Color',basalcolor,'LineWidth',0.1);hold on
% ylim([0 8])
% xlim([0 300])
% subplot(2,1,2);plot(TimeVector,FC(:,frameshift:end),'Color',basalcolor,'LineWidth',0.1);hold on
% ylim([0 8])
% xlim([0 300])


% 
% figure(2)
% subplot(2,2,1);plot(TimeVector,CONCENTRATION(:,frameshift:end),'LineWidth',0.1);
% ylim([0 8])
% xlim([0 120])
% subplot(2,2,3);plot(TimeVector,FC(:,frameshift:end),'LineWidth',0.1);
% ylim([0 8])
% xlim([0 120])
% stophere=1;

figure(2)
subplot(2,1,1);j =plot(TimeVector,FC(:,frameshift:end),'LineWidth',0.1);
ylim([0 8])
xlim([0 120])
subplot(2,1,2);jj =plot(TimeVector,FC(:,frameshift:end),'LineWidth',0.1);hold on
ylim([0 8])
xlim([0 120])
stophere=1;

h{1} = j;
h{2} = jj;
end

function  plotCalledPeaks(PLOT,NAMES,peaktime,x,peakABS,peakDIFF,peakFC,timetopeakFC,tresponseFC,y,fresponseFC)
figure(2)
for i=1:1:length(PLOT);
subplot(3,3,i);plot((0:x:((size(PLOT{i},2)*x)-1)),PLOT{i}); title(NAMES{i});hold on
axis([0 peaktime 0 (max(max(PLOT{i})))]);
%legend('show');
i+1;
end
for j = 1:1:length(peakABS);
subplot(3,3,1);plot((timetopeakFC(j)),peakABS(j),'--rs');hold on
axis([0 peaktime 0 (max(peakABS))*1.2])
subplot(3,3,2);plot((timetopeakFC(j)),peakDIFF(j),'--rs'); hold on
axis([0 peaktime 0 (max(peakDIFF))*1.2])
subplot(3,3,3);plot((timetopeakFC(j)),peakFC(j),'--rs');hold on
subplot(3,3,3);plot(((tresponseFC(j)+(y-1)*x)),fresponseFC(j),'o');hold on
%subplot(3,3,3);plot(((timetopeak50(j))-x),f50peak(j),'o');hold on
%subplot(3,3,3);plot((adaptationtime(j)),ADAPTmin(j),'*');hold on
axis([0 peaktime 0 (max(peakFC))*1.2])
j+1;
end
hold off;
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
end

function [NAMES,PLOT] = timecourseplotone(response,TIME,y,CFP)
%calculate basals
for i = 1:1:size(response,1)
basal(i) = nanmean(response(i,1:1:y));
basalcfp(i) = nanmean(CFP(i,1:1:y));
end

%absolute levels plot
NAMES{1}='f(x)';
PLOT{1} = response;

%absolute difference plot
%calcuate basal vetor for subtraction to generate absolute difference
for i=1:1:size(response,1);
bs = (ones(size(response)));
basalsubtract(i,:)= bs(i,:).*basal(i);
end
absolutedifference = response-basalsubtract;
NAMES{2}= 'f(x)-f(i)';
PLOT{2} = absolutedifference;

%fold change response plot 
%(divide by 2 frames preceeding basal if possible)
for i=1:1:size(response,1);
    if isnan(nanmean([response(i,y),response(i,y-1)]))==1
        FC(i,:) = response(i,:)./basal(i);
    else
        FC(i,:) = response(i,:)./nanmean([response(i,y),response(i,y-1)]);
    end
end
NAMES{3}='f(x)/f(i)';
PLOT{3} = FC;

%smoothed functions plot
%generate smoothed functions
for i= 1:1:size(FC,1)
yyFC(i,:) = smooth(TIME,FC(i,:),0.03,'moving');
end
NAMES{4} = 'smoothFC';
PLOT{4} = yyFC;

%medianFC plot
medianFC = nanmedian(FC);
NAMES{5} = 'medianFC';
PLOT{5} = medianFC;

%median CFP plot
for i = 1:1:size(response,1)
    fccfp(i,:) = CFP(i,:)/basalcfp(i);
end
NAMES{6} ='medianCFP';
PLOT{6} = nanmedian(fccfp) ;

% artifact plot
NAMES{7} = 'artifact';
PLOT{7} = 1;

%gradient plot
for i = 1:size(FC,1)
    respy(i,:) = gradientdefined(FC(i,1:end),TIME(1:end));
end
PLOT{8} = respy;
NAMES{8} = 'FC gradient';


%fold change in CFP plot
for i = 1:1:size(response,1)
    fccfp(i,:) = CFP(i,:)/basalcfp(i);
end
NAMES{9} = 'fcCFP';
PLOT{9} = fccfp;
end

function plotCVs(SCALARS,NAMO,COLORS,CV)
Coefficient = figure(4);

% Create axes
for i = 1:1:length(CV)
    XTickLabel{i}=NAMO{i};
end

CVtouse = [1:length(SCALARS)];
axes1 = axes('Parent',Coefficient,...
 'XTickLabel',XTickLabel(CVtouse),'XTick',[1:1:length(CVtouse)]);    
ylim(axes1,[0 1.1]);
box(axes1,'on');
hold(axes1,'all');
baz=1;
for i = CVtouse
    bar(baz,CV(i),'FaceColor',COLORS{i});hold on
    baz = baz+1;
end
    end

function plotScalars(SCALARS,CV,basal,NAMO)

%% SCALARS PLOTS--SCATTERz
for i=1:1:length(SCALARS);
colourmax(i) = (max(SCALARS{i}/min(SCALARS{i})));
i=i+1;
end
colours = colourmax./max(colourmax);

% figure %plot none correlated single cell data
% for i=1:1:length(SCALARS);
% %subplot(3,3,i);plot((0:1:length(SCALARS{i})),SCALARS{i}); title(NAMO{i});
% %subplot(3,6,i);scatter(1:1:length(SCALARS{i}),SCALARS{i}); title(NAMO{i});
% subplot(4,5,i);scatter(ones(length(SCALARS{i}),1).*1.4,SCALARS{i}); title(NAMO{i});hold on%set(gca,'Color',[0 colours(i) 0]);
% subplot(4,5,i);scatter(1.4,median(SCALARS{i}),'rs');
% %axis([0 length(SCALARS{i}) 0 max(max(SCALARS{i}))]);
% %legend('show');
% i+1;
% end
% hold off;


            CVcolours = CV./max(CV);
            for i=1:1:length(CV);
            if isnan(CVcolours(i)) == 1
            CVcoloursx{i} = [1 0 0];
            else
                CVcoloursx{i} = [0 CVcolours(i) 0];
            end
            i=i+1;
            CVcoloursx;
            end

%% Plot SCALARS vs Basal level--attempt to find correlations
figure(3)
for i=1:1:length(SCALARS);
% subplot(4,5,i);scatter(basal,SCALARS{i});title(NAMO{i});hold on;
subplot(5,6,i);semilogx(basal,SCALARS{i},'LineStyle','none','Marker','o');title(NAMO{i});hold on;
set(gca,'Color',CVcoloursx{i});
i+1;
end
hold off
end



function basalcolor = calcbasalcolor(response,colourz,basal)
for i=1:1:size(response,1);
orange = colourz.orange; blue = colourz.blue; red = colourz.red; yellow = colourz.yellow; green = colourz.green; green = colourz.green; cyan = colourz.cyan; grey = colourz.grey; black = colourz.black;

medbasal=median(basal);
basalm.m200 = medbasal*2.00;
basalm.m175= medbasal*1.75;
basalm.m150 = medbasal*1.5;
basalm.m125 = medbasal*1.25;
basalm.m100 = medbasal;
basalm.m075 = medbasal*0.75;
basalm.m050 = medbasal*0.50;
basalm.m025 = medbasal*0.25;



%% basal colors
if basal(i) > basalm.m200
    basalcolor{i} = grey;
elseif basal(i) > basalm.m175
basalcolor{i} = black;
elseif basal(i) > basalm.m150
basalcolor{i} = green;
elseif basal(i) > basalm.m125
 basalcolor{i} = red;
elseif basal(i) > basalm.m100;
    basalcolor{i} = orange;
elseif basal(i) > basalm.m075
    basalcolor{i} = yellow;
elseif basal(i) > basalm.m050
    basalcolor{i} = blue;
elseif basal(i) > basalm.m025
    basalcolor{i} = green;
else
    basalcolor{i} = cyan;
end
end
end

function plotthatshadederrorbar(FC,medcorr,x)      
%plotshadederrorbarstuff
    meanFC = nanmean(FC);
    medianFC = nanmedian(FC);
    stdFC = nanstd(FC);
    NAMES{4} = 'meanFC';
    PLOT{4} = meanFC;
    FCerrorup = meanFC+(stdFC/2);
    FCerrordown = meanFC-(stdFC/2);

if size(meanFC,2)>1
figure(5)
subplot(2,1,1);shadedErrorBar(0:x:(size(medcorr,2)-1)*x,medcorr,{@nanmedian,@nanstd},'b',1);hold on
subplot(2,1,2);shadedErrorBar(0:x:(size(meanFC,2)-1)*x,FC,{@nanmedian,@nanstd},'b',1);hold on
hold off
end
end

function graddefined= gradientdefined(vector,time)
for j = 1:size(vector,2)
    if j == 1
risedefined(1,j) = (vector(1,j+1)-vector(1,j));
rundefined(1,j) = (time(1,j+1)-time(1,j));
    elseif j == size(vector,2)
risedefined(1,j) = (vector(1,j)-vector(1,j-1));
rundefined(1,j) = (time(1,j)-time(1,j-1));      
    else
risedefined(1,j) = (0.5).*(vector(1,j+1)-vector(1,j-1));
rundefined(1,j) = (0.5).*(time(1,j+1)-time(1,j-1));
    end
graddefined(1,j) = risedefined(1,j)./rundefined(1,j);
end
end

function [peakVAR,tpeakVAR,timetopeakVAR] = callpeak(matrix,peaklength,y,x)
for i = 1:size(matrix,1)
peakz(i) = max(matrix(i,y:1:peaklength));
tpeakVAR(i,:) = find(matrix(i,1:1:peaklength)==peakz(i),1,'first');
timetopeakVAR(i) = (tpeakVAR(i,:)-1)*x;
peakVAR(i,:) = nanmean([matrix(i,tpeakVAR(i,:)),matrix(i,tpeakVAR(i,:)+1),matrix(i,tpeakVAR(i,:)+2)]);
end
end

function [gradVAR,rateVAR,tresponseVAR,fresponseVAR] = computegradient(smoothedmatrix,TIME,peaklength,y,x,tshift)
for i = 1:size(smoothedmatrix,1)
respy(i,:) = gradientdefined(smoothedmatrix(i,1:end),TIME(1:end));
resp(i,:) = gradientdefined(smoothedmatrix(i,y:peaklength),TIME(y:peaklength));
varrespo(i) = max(resp(i,2:end-1));
if isnan(varrespo(i))
varrespon(i) = NaN;
gradVAR(i,:) = NaN;
rateVAR(i) = NaN;
    
    tresponVAR(i) = NaN;
%     tresponseVAR(i) = tresponVAR(i)*x;
    tresponseVAR(i) = NaN;
    fresponseVAR(i) = NaN;
    disp('lostone')
else
varrespon(i) = find(resp(i,:) == varrespo(i),1,'first');
gradVAR(i,:) = gradientdefined(smoothedmatrix(i,y:peaklength),TIME(y:peaklength));
rateVAR(i) = mean([gradVAR(i,varrespon(i)),gradVAR(i,varrespon(i)+1),gradVAR(i,varrespon(i)-1)]);

    
    tresponVAR(i) = find(resp(i,:) == varrespo(i),1,'first');
%     tresponseVAR(i) = tresponVAR(i)*x;
    tresponseVAR(i) = TIME(y+tresponVAR(i)-1)-(TIME(y)+tshift);
    fresponseVAR(i) = smoothedmatrix(i,y+tresponVAR(i));
end
end
end


function [SCALARS,NAMO,COLORS] = timeResponses(response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
t = time;
% for sn = sno:sno+14;
% SCALARS{sn} = response(:,y+t);
% NAMO{sn} = strcat(responseName,num2str(t*x));
% COLORS{sn} = color;
% t = t+1;
% end
% for sn = sno:sno+14;
tlength = (peakLEngth./x)-time;
for sn = sno:sno+tlength;
        SCALARS{sn} = response(:,y+t);
NAMO{sn} = strcat(responseName,num2str(t*x));
COLORS{sn} = color;
t = t+1;
end

end

function [SCALARS,NAMO,COLORS] = timeResponsesAveraged(response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
t = time;
% for sn = sno:sno+14;
% SCALARS{sn} = response(:,y+t);
% NAMO{sn} = strcat(responseName,num2str(t*x));
% COLORS{sn} = color;
% t = t+1;
% end
LOAvg = 3;
tlength = (peakLEngth./x)-time;
for sn = sno:sno+tlength;
    if y+t+5> size(response,2)
        SCALARS{sn} = nanmean(response(:,y+t:end),2);
    else
        SCALARS{sn} = nanmean(response(:,y+t:y+t+LOAvg),2);
    end
NAMO{sn} = strcat(responseName,num2str(t*x),'ThreeFrameAvg');
COLORS{sn} = color;
t = t+1;
end

end

function [SCALARS,NAMO,COLORS] = timeResponsesRate(response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
t = time;

% tlength = size(response,2)-(y+3)-time;
tlength = (peakLEngth./x)-time-2;
for sn = sno:sno+tlength;
% for sn = sno:sno+12;
SCALARS{sn} = response(:,y+t);
NAMO{sn} = strcat(responseName,num2str(t*x));
COLORS{sn} = color;
t = t+1;
end
end


function [SCALARS,NAMO,COLORS] = timeResponsesIntegrated(response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
t = time+2;
% tlength = size(response,2)-(y+3)-time;
tlength = (peakLEngth./x)-time-2;
for sn = sno:sno+tlength;
if y+t>size(response,2)
Hs = nanmean(response(:,y:y+1),2);
Hn = Hs.*NaN;
SCALARS{sn} =Hn;
else
SCALARS{sn} = nanmean(response(:,y:y+t),2);
end
NAMO{sn} = strcat(responseName,num2str(t*x));
COLORS{sn} = color;
t = t+1;
end
end

function [SCALARS,NAMO,COLORS] = timeResponsesEndAveraged(TIME,response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
t = time;

% for i= 1:1:size(response,1)
% yyResponse(i,:) = smooth(TIME,response(i,:),0.03,'moving');
% end

tlength = (peakLEngth./x)-time-4;
for sn = sno:sno+tlength;
% if y+15>size(response,2)
SCALARS{sn} = nanmean(response(:,y+t:end),2);
% else
% SCALARS{sn} = nanmean(response(:,y+t:y+t+4),2);
% end
NAMO{sn} = strcat(responseName,num2str(t*x),'avgToEnd');
COLORS{sn} = color;
t=t+1;
end

% figure(33);plot(TIME,response');hold on


end


function [SCALARS,NAMO,COLORS] = maxVarOvertime(TIME,response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
t = time+0;

% for i= 1:1:size(response,1)
% yyResponse(i,:) = smooth(TIME,response(i,:),0.03,'moving');
% end

% tlength = size(response,2)-(y+3)-time;
tlength = (peakLEngth./x)-time-2;
for sn = sno:sno+tlength
    for i = 1:size(response,1)
maxx(i) = nanmax(response(i,y-1:y+t));
    end
SCALARS{sn} = maxx;
NAMO{sn} = strcat(responseName,num2str(t*x),'t');
COLORS{sn} = color;
t=t+1;
end
% figure(33);plot(TIME,response');hold on


end




function [SCALARS,NAMO,COLORS] = timeFifty(response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
sn=sno;
t=time+0;
tlength = (peakLEngth./x)-time-2;
responsetrunc = response(:,y+1:end); %find t50 for whole trace
responsetrunc = response(:,y+1:y+tlength); %find t50 for just the "first peak"
maxx = max(responsetrunc,[],2);
for i = 1:size(responsetrunc,1)
    idx = find(responsetrunc(i,:) > maxx(i)./2,1,'first');
    if isempty(idx)
       idx = NaN;
    end
    
    time(i) = (idx-1).*x;
end    
    


    SCALARS{sn} = time;
    NAMO{sn} = strcat(responseName);
    COLORS{sn} = color;
    t=t+1;
    % figure(33);plot(TIME,response');hold on
    
end



function [SCALARS,NAMO,COLORS] = ResponsePostTime(response,SCALARS,NAMO,COLORS,sno,time,responseName,color,x,y,peakLEngth)
sn=sno;
t=time-4;

responsetime = SCALARS{2};
rt = round((responsetime./x)+1);

tlength = (peakLEngth./x)-time-0;
for sn = sno:sno+tlength
    
    for i = 1:size(response,1)
        if y+rt(i)+t>size(response,2)
        fcpostRT(i) = response(i,end);    
        elseif y+rt(i)+t<y
            fcpostRT(i) = response(i,y);
        else
        fcpostRT(i) = response(i,y+rt(i)+t);
        end
    end

    SCALARS{sn} = fcpostRT;
    NAMO{sn} = strcat(responseName,num2str(t.*x));
    COLORS{sn} = color;
    t=t+1;
end

end