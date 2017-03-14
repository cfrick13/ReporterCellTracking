
%*files must be exported as TIFFs from Zen using channel names
%*the export folder must have the same name as the parent folder containing the zen experiment files
% such as parent folder = 'D:/kkim/2016_10_25 bcat cells' then export folder should be 'D:/kkim/2016_10_25 bcat cells/2016_10_25 bcat cells'

% parentFolderName = 'FrickPaperData';
% parentFolderName = 'kkim';
% mfile = mfilename('fullpath');
% [a,b] = regexp(mfile,parentFolderName);
% mfiledir = mfile(1:b+1);
% parentdir = mfiledir;
% A = strcat(parentdir);


%set parent directory
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
% parentdir = mdir(1:b(end-1));
% exportdir = strcat(parentdir,'TrackingMstack/Export/');
parentdir = mdir(1:b(end));
exportdir = strcat(parentdir,'Export/');





%for loop iterating through subdirectories
% for BB = {'2017_01_22 plate','2017_01_25 plate','2017_01_27 plate'};
% for BB = {'2017_01_30 plate exp2','2017_01_27 plate exp1','2017_01_31 plate exp1'};
% for BB = {'2017_01_27 plate exp1','2017_01_31 plate exp1'};
cd(parentdir)

% for BB = {'2014_09_19 plate exp1','2014_09_27 plate exp1','2014_09_30 plate exp1','2014_10_01 plate exp1','2014_10_04 plate exp1','2014_11_21 plate exp1'}
% for BB = {'2014_09_27 plate exp1','2014_09_30 plate exp1','2014_10_01 plate exp1','2014_10_04 plate exp1','2014_11_21 plate exp1'}    
for BB = {'2017_03_13 plate exp1'}

% for BB = {'2017_03_02 plate exp1'}
B = char(BB);
FileName = B;
disp(FileName)

supertic = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExtractMetadataAndImages(B);
% ExtractMetadataOnly(B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datequery = strcat(FileName,'*DoseAndScene*');
cd(exportdir)
filelist = dir(datequery);
    if isempty(filelist)
        error(strcat('need to run ExtractMetadata for-',FileName));
%        dosestruct = makeDoseStruct; %run function to make doseStruct 
    else
        dosestructstruct = load(char(filelist.name));
        dosestruct = dosestructstruct.dosestruct;
    end
A=[];
channelstoinput = dosestructstruct.channelNameSwapArray;
bkg = dosestructstruct.BACKGROUND;
BACKGROUND = bkg{1};
dontsegment = BACKGROUND;
% dontsegment = horzcat(1:5,BACKGROUND);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);
totalTime = toc(supertic);
disp(strcat('total time from extract to flat is=', num2str(totalTime./60),' minutes'));
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SortBKGandFLATcorrectedImagesIntoFolders(A,B,channelstoinput);
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BleachCorrectionOfTimeLapseImages(A,B,channelstoinput);
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SegmentationOfTimeLapseImages(A,B,dontsegment);
% renameDirectoryItems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end



% A = 'D:\Frick\';
% B = '2016_01_25 smad3g smFISH';
% % 
% % preprocessBKGandFLATnew(A,B)
% % % lineageyo(A,B);
% % lineageyoz(A,B);
% % cd('D:\Users\zeiss\Documents\MATLAB')
% % % bleachcorrection(A,B)
% % bleachcorrectionNew(A,B)
% cd('C:\Users\zeiss\Documents\MATLAB')
% toughsegmentationforstacks(A,B)

