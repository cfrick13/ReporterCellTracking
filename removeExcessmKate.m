function removeExcessmKate(A,B)
%%
% channelinputs = '(mKate|EGFP|CFP|DIC)';
channelstoinput = {'mKate','EGFP','CFP','DIC'};
% channelstoinput = {'_mKate','_EGFP','_CFP','DIC'};

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


cd('mKate_flat')
filelist = dir('*CFP*')
if ~isempty(filelist)
delete(filelist.name);
end

filelist = dir('*EGFP*')
if ~isempty(filelist)
delete(filelist.name);
end

filelist = dir('*DIC*')
if ~isempty(filelist)
delete(filelist.name);
end
stophere=1;

cd ..
cd ..
end