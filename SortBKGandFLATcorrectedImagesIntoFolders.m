function SortBKGandFLATcorrectedImagesIntoFolders(A,B)
%%
% channelinputs = '(mKate|EGFP|CFP|DIC)';
channelstoinput = {'mKate','EGFP','CFP','DIC'};
% channelstoinput = {'_mKate','_EGFP','_CFP','DIC'};

channelinputs = '(';  % result is channelinputs = '(mKate|EGFP|CFP|DIC)';
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
    for folder = {folderlist.name} %set the scene directory for all
    cd(folder{1})
    mkdir('tiffs');
    %check if directories have been made previously
    %make individual directories
        for i = 1:length(channelstoinput)
            dd = dir(strcat(channelstoinput{i},'_flat'));
                if isempty(dd) 
                    nodir = 1;
                    filelist = dir(strcat('*_',channelstoinput{i},'*.tif')); 
                        if isempty(filelist) %if no files are present do not make the directory
                        else
                        mkdir(strcat(channelstoinput{i},'_flat'))
                        end
                else
                    nodir = 0;
                end
        
            foldername = strcat(channelstoinput{i},'_flat');
            filelist = dir(strcat('*_',channelstoinput{i},'*.tif'));  
                if isempty(filelist) && nodir==0
                    cd(foldername)
                    filelist = dir(strcat('*_',channelstoinput{i},'*.tif'));  

                    if isempty(filelist) 
                        cd .. 
%                         disp(foldername)
                        ddd = dir(foldername);
                        if isempty(ddd)
                        else
                            disp(foldername)
                            disp(ddd)
                            pause(0.5)
                        rmdir(foldername);
                        end
                    else
                        cd ..
                    end

                else

                    for cfile={filelist.name}
                       filepath = which(cfile{1});
                       disp(filepath);
                       movefile(filepath,foldername)
                    end
                end
        end
        cd ..
    end
end
