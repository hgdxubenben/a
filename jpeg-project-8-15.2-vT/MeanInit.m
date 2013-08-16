function meanI = MeanInit(filename, CoveredBuf)
%?算已恢复文件部分的ED值
outfile = fopen(filename,'wb');  
while outfile == -1        
    outfile = fopen(filename,'wb');        
end

fwrite(outfile,CoveredBuf);
fclose(outfile);

rgb = imread(filename); %the unompressed data
Location = FindLocationRecovered(filename);  %the impressed data 
row = Location.row;
col = Location.col;
width = size(rgb, 2);
meanI = 0;
for i = 2: row
    upR   = rgb(i-1, :,:);
    midR  = rgb(  i, :,:);
    value = fdiffMU_1D(midR, upR)/width;
    meanI = meanI + value;
end;
meanI; %the average similarity between lines
meanI = meanI / (row-1);
