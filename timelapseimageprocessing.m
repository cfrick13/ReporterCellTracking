% A = 'D:\Users\zeiss\Pictures\Frick\';
A = 'D:\Users\zeiss\Pictures\KKIM\';
% B = '2016_05_31 screen';
B = '2016_06_08 bcat-cit';
% renamemCherrytoMkate(A,B)
cd('D:\Users\zeiss\Documents\MATLAB')
preprocessBKGandFLATofTimeLapseImages(A,B);
cd('D:\Users\zeiss\Documents\MATLAB')
SortBKGandFLATcorrectedImagesIntoFolders(A,B);
cd('D:\Users\zeiss\Documents\MATLAB')
BleachCorrectionOfTimeLapseImages(A,B);
cd('D:\Users\zeiss\Documents\MATLAB')
% SegmentationOfTimeLapseImages(A,B);
SegmentationOfTimeLapseImagesHoechst(A,B);

%run uiTrackThatCells
