function lineageyozX(A,B)
%%
F = {'_c1_ORG.tif','_c2_ORG.tif','_c3_ORG.tif'};
FF = {'*c1_*.tif','*c2_*.tif','*c3_*.tif'};
cd (strcat(A,B));

cd ('flatfield_corrected')

%%
folderlist = dir('*_s*');
folderlist.name;

for folder = {folderlist.name}
cd(folder{1})
mkdir('c1_flat')
mkdir('c2_flat')
mkdir('c3_flat')
mkdir('tiffs');

filelist = dir('*c2*.tif');  
for cfile={filelist.name}
       filepath = which(cfile{1});
%        [a,b,~,~] = regexp(filepath,'c1');
%        c = filepath(a:b);
%        c(2) = '2';
%        filepath(a:b) = c;
       disp(filepath)
       movefile(filepath,'c1_flat')
end

filelist = dir('*c1*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
%        [a,b,~,~] = regexp(filepath,'c2');
%        c = filepath(a:b);
%        c(2) = '1';
%        filepath(a:b) = c;
       disp(filepath);
       movefile(filepath,'c2_flat')
end

filelist = dir('*c3*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'c3_flat')
end
 cd .. 

end
end