function SortBKGandFLATcorrectedImagesIntoFolders(A,B)
%%
% channelinputs = '(mKate|EGFP|CFP|DIC)';
channelstoinput = {'mKate','EGFP','CFP','DIC'};
channelinputs = '(';
for i=1:length(channelstoinput)
    if i ==1
    channelinputs = strcat(channelinputs,channelstoinput{i});
    elseif i < length(channelstoinput)
        channelinputs = strcat(channelinputs,'|',channelstoinput{i});
    else
        channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
    end
end

cd (strcat(A,B));
cd ('flatfield_corrected')

%%
folderlist = dir('*_s*');
folderlist.name;

for folder = {folderlist.name}
cd(folder{1})

for i = 1:length(channelstoinput)
mkdir(strcat(channelstoinput{i},'_flat'))
end
mkdir('tiffs');

filelist = dir('*mKate*.tif');  
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'mKate_flat')
end

filelist = dir('*EGFP*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'EGFP_flat')
end

filelist = dir('*CFP*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'CFP_flat')
end

filelist = dir('*DIC*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'DIC_flat')
end
 cd .. 

end
end