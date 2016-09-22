function renamemCherrytoMkate(A,B)
C = B;

wrongname = {'*HEX*'

E={B}; %beginning of file name for saving
channelstoinput = {'_mKate','_EGFP','_CFP','DIC'};
% channelstoinput = {'mKate','_EGFP','_CFP','_DIC'};
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



BACKdir = strcat(A,B,'\back');mkdir (BACKdir);
DIRone = strcat(A,B,'\',C);

cd(DIRone);
dirlist = dir('*mCherry*.tif');
cfile = {dirlist.name};
parfor i = 1:length(cfile)
    oldfilename = cfile(i);
    oldfile = oldfilename{1};
    [a,b] = regexp(oldfile,'mCherry');
    newfile = strcat(oldfile(1:a-1),'mKate',(oldfile(b+1:end)));
    movefile(oldfile,newfile);
    disp(i)
end


end