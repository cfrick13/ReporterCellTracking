% A = 'D:\Frick\';
% for BB = {'2014_03_15 smad3g smFISH'};
% for BB = {'2016_05_31 screen','2016_06_01 screen'};
for BB = {'2016_10_12 plate exp4 screen'};
% for BB = {'2016_10_03 clone+32 part1','2016_10_03 clone+32 part2'};
% for BB = {'2015_08_31 smad3g smFISH','2015_09_03 smad3g smFISH','2015_12_15 smad3g smFISH','2015_12_19 smad3g smFISH','2016_01_25 smad3g smFISH','2016_02_09 smad3g smFISH','2016_02_19 smad3g smFISH'};
B = BB{1};
% A = 'D:\Frick\2016_09_28 screen\';
% A = 'D:\Frick\2016_10_01 plate -- clone+32 and hoechst imaging\';
A = 'D:\Frick\2016_10_12 plate\';
% A = 'G:\';
% A = 'D:\Users\zeiss\Pictures\Frick\';
% B = '2016_05_31 screen';
% B = '2016_06_04 screen';
% if strcmp(B,'2016_06_01 screen')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% renamemCherrytoMkate(A,B)
% renamemWRONGtoRIGHT(A,B)
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preprocessBKGandFLATofTimeLapseImages(A,B);
cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SortBKGandFLATcorrectedImagesIntoFolders(A,B);
cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BleachCorrectionOfTimeLapseImages(A,B);
cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SegmentationOfTimeLapseImages(A,B);
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

