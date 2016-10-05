function AllCellTimeTracesMakerMacWin
ExperimentStructure = struct();

% for Dates = {'2014_09_30','2014_10_01','2014_10_04','2014_10_24','2014_11_25','2015_01_15','2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01','2015_08_31','2015_09_03'};  
% for Dates = {'2014_09_30','2014_10_01','2014_10_04','2014_10_24','2014_11_25','2015_01_15','2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01','2015_07_02','2015_07_10','2015_07_14','2015_08_31','2015_09_03'};    
% for Dates = {'2014_09_30','2014_10_01','2014_10_04','2014_10_24','2014_11_25','2015_01_15','2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01'};
% for Dates = {'2015_08_31','2015_09_03'};
% for Dates = {'2014_09_30','2014_10_01','2014_10_04'}
    % for Dates = {'2015_07_02','2015_07_10','2015_07_14'}
    % for Dates = {'2015_07_02'};    
% for Dates = {'2015_08_31','2015_09_03'};  
% for Dates = {'2015_12_15','2015_12_19'}
for Dates = {'2016_10_01'}
Date = char(Dates);
ADir = 'D:\Users\zeiss\Pictures\Frick\ImagingResults';
BDir = 'D:\Users\zeiss\Documents\';
cd(ADir)
cd(ADir)
DetailsFilename = 'ExperimentalDetails.xlsx';
% DetailsFilename = 'ExperimentalDetailsNoTimeAdjust.xlsx';
[detailssheet,detailstext,detailsraw] = xlsread(DetailsFilename,Date);

filename = strcat(Date,'*smad3g same cells*SHEET.xlsx');
filelist = dir(filename);
[status,sheets] = xlsfinfo(char({filelist.name}));

DateStruct = struct();
Neato = struct();
for sheet = sheets
    cd(ADir)
    svalue = char(sheet);
    [aa,bb] = regexp(svalue,'s[0-9]+');
    svaluestring = svalue(aa:bb);
    if isempty(svaluestring)
    else
        
        
        %     Date
if strcmp(Date,'2015_03_31')
    if strcmp(svalue,'s27')
    stophere=1;
    end
end       
        
    deta = createTimeTraceFile(svalue,filelist,detailssheet,detailsraw);
    

    
    
    cd(strcat(BDir,'MATLAB\AllImagingDataCompiled\Newz'))
    Dose = char(deta.dose);
    Frame = deta.frame;
    frame = num2str(Frame);
    if size(frame,2)==1
    Frame = '00';
    Frame(2) = frame;
    else
    Frame = frame;
    end
    pvalue = svalue;
    pvalue(1) = 'p';
    
    if strcmp(pvalue,'p18')
        stophere=1;
    end
    DATA = vertcat(deta.timeVector,deta.cellTraces);
    savename = strcat('data_',Date,'_smad3g_dose',Dose,'_frame',Frame,'_',pvalue)
    Neato.(savename) = DATA;
    save(strcat(savename,'.mat'),'-struct','Neato',savename);
%     save(strcat(char(emmiesname),'.mat'), '-struct','M',emmiesname);
%     data_2014_09_30_smad3g_dose0dot00_frame08_p25
    
    
    DateStruct.(svalue) = deta;
    end
end
ExperimentStructure.(strcat('i',Date,'i')) = DateStruct;

end

cd(strcat(BDir,'MATLAB\AllImagingDataCompiled'))

save('ExperimentStructure.mat','ExperimentStructure');
% save('ExperimentStructureNoTimeAdjustMat.mat','ExperimentStructure');
stophere=1;
end


function details = createTimeTraceFile(svalue,filelist,detailssheet,detailsraw)

[spreadsheet,spreadtext,spreadraw] = xlsread(char({filelist.name}),svalue);

% nanidx = cellfun(@(x) sum(isnan(x)),spreadraw,'UniformOutput',1);
% CC = bwconncomp(~nanidx);
% PX = CC.PixelIdxList;
% imgsize = CC.ImageSize;
%   
% %minrow
% %mincol
%     for j=1:length(PX)
%         px = PX{j};
%         [row,col] = ind2sub(imgsize,px);
%         minrow(j) = min(row);
%         maxrow(j) = max(row);
%         mincol(j) = min(col);
%         maxcol(j) = max(col);
%     end
% 
%     
% %find cell comments
% numidx = cellfun(@(x) isnumeric(x),spreadraw,'UniformOutput',1);
% CC = bwconncomp(~numidx);
% PX = CC.PixelIdxList;
% imgsize = CC.ImageSize;
%  
% %minrownum
% %mincolnum
%     for j=1:length(PX)
%         px = PX{j};
%         [row,col] = ind2sub(imgsize,px);
%         minrownum(j) = min(row);
%         maxrownum(j) = max(row);
%         mincolnum(j) = min(col);
%         maxcolnum(j) = max(col);
%     end
% 
% if length(PX)<1
%     nocomments=1;
%     commentcolnum=1;
% else
%     nocomments=0;
%     commentcolnum=mincolnum(end);
% end
%     
%     
% %find the cell traces values (background and cells are the two pixel areas) 
% % nanidx = 
% idxforvalues = ~nanidx & numidx;
% CC = bwconncomp(idxforvalues);
% PX = CC.PixelIdxList;
% imgsize = CC.ImageSize;
%   
%     cycle=1;
%     for j=1:length(PX)
%         px = PX{j};
%         if length(px)>(maxrow(end))
%             PXX{cycle} = px;
%             cycle=cycle+1;
%             disp(length(px))
%         end
%     end
%     CC.PixelIdxList = PXX;
%     CC.NumObjects = length(PXX);
%     
%     PX = CC.PixelIdxList;
%     imgsize = CC.ImageSize;
%     for j=1:length(PX)
%         px = PX{j};
%         [row,col] = ind2sub(imgsize,px);
%         minrowval(j) = min(row);
%         maxrowval(j) = max(row);
%         mincolval(j) = min(col);
%         maxcolval(j) = max(col);
%     end
%     
%     removeidx = (mincolval==maxcolval);
%     
%         if sum(removeidx)>0
%             removeidxnum = find(removeidx==1);
%             mincolval(removeidxnum)=[];
%             maxcolval(removeidxnum)=[];
%             maxrowval(removeidxnum)=[];
%             minrowval(removeidxnum)=[];
%         end
%         
%         
% if numel(maxcolval)<2
%     stophere=1;
% end
% if diff([maxcolval(1) mincolval(1)])>diff([maxcolval(2) mincolval(2)])
%     img = ~nanidx;
%     se = [0 1 0; 0 1 0; 0 1 0];
%     stophere=1;
%     imagesc(img);
%     imgc = imopen(img,se);
%     imagesc(imgc);
% end
%         
% backAverages(1,:) = nanmean(spreadsheet(minrowval(1):2:maxrowval(1),mincolval(1):maxcolval(1)),1);
% backAverages(2,:) = nanmean(spreadsheet(minrowval(1)+1:2:maxrowval(1),mincolval(1):maxcolval(1)),1);
%     if sum(size(spreadsheet) == size(spreadraw))==2
%         stophere=1;
%     else
%         stophere=1;
%         disp('poop')
%     end
%     if numel(maxrowval)==1
%        stophere=1; 
%     end
%     
%     if strcmp(svalue,'s19')
%         stophere=1;
%     end
% rawTraces = zeros(maxrowval(2)-minrowval(2)+1,maxcolval(2)-mincolval(2)+1);
% rawTraces(1:2:(maxrowval(2)-minrowval(2)),:) = spreadsheet(minrowval(2):2:maxrowval(2),mincolval(2):maxcolval(2));
% rawTraces(2:2:(maxrowval(2)-minrowval(2)+1),:) = spreadsheet(minrowval(2)+1:2:maxrowval(2),mincolval(2):maxcolval(2));
% 
% cellTraces = vertcat(backAverages,rawTraces);
% if nocomments==1
%     cellComments = cell(maxrowval(end)+1-minrowval(end),1);
% else
% cellComments = spreadraw(minrowval(end):maxrowval(end),commentcolnum);
% end
% 
% cellTracesLength = size(cellTraces,2); 
% cellNumbers = spreadsheet(minrow(2):2:maxrow(2),mincol(2));



spreadBacks = spreadsheet(3:14,:);
kkk = find(spreadBacks<10 & spreadBacks>0);
[yy,zz] = ind2sub(size(spreadBacks),kkk);

if isempty(zz)
    endhere=1;
end

backnumbers = spreadBacks(yy,zz(1));
backAverages = spreadBacks(yy(end)+1:yy(end)+2,zz(1)+2:end);

spreadTraces = spreadsheet(10:end,:);
jjj = find(spreadTraces<10 & spreadTraces>0);
[a,b] = ind2sub(size(spreadTraces),jjj);
cellnumbers = spreadTraces(a(1):end,b(1));

sizedifference = size(spreadraw,2)-size(spreadTraces,2);

spreadraw(a(1)+9-1,b(1)+sizedifference) = {'CELLS'};
spreadraw(a(1)+9-1,b(1)+1+sizedifference) = {'COMMENTS'};


disp(svalue)
if strcmp(svalue,'s18')
    stophere=1;
end
if (zz(1)-b(1))<0
rawTraces = spreadsheet(a(1)+9:end,b(1)+2:end);
else
rawTraces = spreadsheet(a(1)+9:end,b(1)+2:end-(zz(1)-b(1)));
end
cellNumbers = spreadsheet(a(1)+9:2:end,b(1));



if size(backAverages,2) == size(rawTraces,2)
cellTraces = vertcat(backAverages,rawTraces);
else
    resizer = backAverages(~isnan(backAverages));
    backAverages = reshape(resizer,[2 size(rawTraces,2)]);
    cellTraces = vertcat(backAverages,rawTraces);
end

cellTracesLength = size(cellTraces,2);  
COMMENTloc = strcmp(spreadraw,'COMMENTS');
[c,d] = ind2sub(size(spreadraw),find(COMMENTloc ==1));
cellComments = spreadraw(c+1:end,d);

[timeVector,details] = timeVectorFunction(cellTracesLength,svalue,filelist,detailssheet,detailsraw);

details.cellNumbers = cellNumbers;
details.cellComments = cellComments;
details.timeVector = timeVector;
details.cellTraces = cellTraces;
details.positon = svalue;
end
    


function [timeVector,details] = timeVectorFunction(cellTracesLength,svalue,filelist,detailssheet,detailsraw)
svalueloc = strcmp(detailsraw,svalue);
[a,b] = ind2sub(size(detailsraw),find(svalueloc ==1));

timeStartloc = strcmp(detailsraw,'TIMESTART');
[c,d] = ind2sub(size(detailsraw),find(timeStartloc ==1));
%use f to determine the column location
%f=2
frameloc = strcmp(detailsraw,'FRAME');
[e,f] = ind2sub(size(detailsraw),find(frameloc ==1));
doseloc = strcmp(detailsraw,'DOSE');
[g,h] = ind2sub(size(detailsraw),find(doseloc ==1));
wellNumberloc = strcmp(detailsraw,'WELL');
[i,j] = ind2sub(size(detailsraw),find(wellNumberloc==1));
fixTimeloc = strcmp(detailsraw,'FIXTIME');
[k,l] = ind2sub(size(detailsraw),find(fixTimeloc==1));
timeAdjustloc = strcmp(detailsraw,'TIMEADJUST');
[m,n] = ind2sub(size(detailsraw),find(timeAdjustloc==1));

%use b to determine the row location corresponding to a specific svalue
%SCENE  + 	TIMESTART	FRAME	DOSE	WELL	FIXTIME	TIMEADJUST
   %b   +     %c          %e     %g       %i       %k       %m

timeInterval = 4;

timeStart   = cell2mat(detailsraw(a,d));
frame       = cell2mat(detailsraw(a,f));
dose        = cellstr(detailsraw(a,h));
well        = strcat('well',num2str(cell2mat(detailsraw(a,j))));
fixTime     = cell2mat(detailsraw(a,l));
timeAdjust  = cell2mat(detailsraw(a,n));

details.timeStart = timeStart;
details.frame = frame;
details.dose = dose;
details.well = well;
details.fixTime = fixTime; 
details.timeAdjust = timeAdjust;



timeVector = zeros(1,cellTracesLength);
timeVector(1:frame) = timeStart:timeInterval:(((frame-1).*timeInterval)+timeStart);
% timeVector(1:frame) = timeStart:timeInterval:(((frame).*timeInterval)+timeStart);
    mid = (frame*timeInterval)+timeAdjust;
    endpoint = ((cellTracesLength-1).*timeInterval)+timeAdjust;
timeVector(frame+1:cellTracesLength) = mid:timeInterval:endpoint;
%need an array with pvalues
%need a matrix that has initial time values in slots that correspond to index positions of pvalues in array


details.timeStart = timeStart;
details.frame = frame;
details.dose = dose;
details.well = well;
details.fixTime = fixTime; 
details.timeAdjust = timeAdjust;


end