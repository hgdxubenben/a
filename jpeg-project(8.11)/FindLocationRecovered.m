function Lo=FindLocationRecovered(filename)

% -- find the end pixel location
global timefind;
global Lotimes;
global sample_factor_v;
global sample_factor_h;
Lotimes = Lotimes + 1;
btime = cputime;
outfile = fopen(filename);  
while outfile == -1         
    outfile = fopen(filename);        
end



sub_fin = imread(filename);     %sub_file_in
width = size(sub_fin, 2);
height = size(sub_fin, 1);

if (mod(width,sample_factor_v*8) ~= 0)      %calculate this when: picture's resolution is not times of 8
    width = width - mod(width, sample_factor_v*8);
    sub_fin = sub_fin(:,1:width,:);
end;
if (mod(height,sample_factor_h*8) ~= 0)
    height = height - mod(height, sample_factor_h*8);
    sub_fin = sub_fin(1:height,:,:);
end;

% sample_factor_v = 2;
% sample_factor_h = 2;

%灰色的RGB=(128,128,128),以此?特征找到?像的?止位置
row_no = flipud(sub_fin(:,end,1));  %提取出?像的最后一列的R值，?列倒置
row    = find(row_no ~= 128, 1);
row    = height - row + 1;
if mod(row,(sample_factor_h*8))~= 0 
    row = round(row/(sample_factor_h*8)) * (sample_factor_h*8);
end;

if row == height %到?文件底部最后的MCU行
    column = width;
else
    col_no = fliplr(sub_fin(row+3,:,1));
    column = find(col_no ~= 128, 1);
    if sum(column) == 0 
        disp('sum column = 0');    
        column = width;
    end
% %     fprintf('row=%d, col=%d\n', row, column);
% %     break;
    column = width - column + 1;
    column = mod((column - sample_factor_v*8 + width),width);
end
column = column - sample_factor_v*8;   %?理?有被完全解析完的MCU
if mod(column,sample_factor_v*8)~= 0 
    column = (floor(column/(sample_factor_v*8)))*(sample_factor_v*8);
end;
if column <= 0
    column = width + column;
end;
if column >= (width - 2*sample_factor_v*8)
    row = row - sample_factor_h*8;
end;

Lo.row = row;
Lo.col = column;
stream = sub_fin( row-(sample_factor_h)*8+1:row, column+1: width,:);
if column ~= width
    stream = cat(2, stream, sub_fin(row+1 : row+(sample_factor_h)*8, 1: column,:));
else
    stream = sub_fin(row+1 : row+(sample_factor_h)*8,1:width,:);
end
Lo.stream = stream;
Lo.len = size(stream,2);
etime = cputime - btime;
timefind = timefind + etime;
% Lo.file = sub_fin;
% Lo
% % % 
% % % row = row + 1;
% % % column = column + 1;
% % % 
% % % Length = floor((row - Beginrow)/8 - 1)*width + (width-Begincol-1) + column;
% % % Lo.len = Length;
% % % 
% % % 
% % % 
% % % % ---- combine and find the data stream
% % % % * * * * combine
% % %     if Lo.row ~= Beginrow
% % %         stream = sub_fin( Beginrow+1 : Beginrow+8, Begincol+1 : width,:);
% % %         for i= (Beginrow/8 +2): (row/8)
% % %             tmp = sub_fin((i-1)*8+1 : i*8,:,:);
% % %             stream = cat(2,stream, tmp);
% % %         end
% % %         tmp = sub_fin(Lo.row+1: Lo.row+8, 1:Lo.col, :);
% % %         stream = cat(2, stream, tmp);
% % %         Lo.stream = stream;
% % %     else
% % %         stream = sub_fin( Beginrow+1 : Beginrow + 8, Begincol+1 : Lo.col,:);
% % %         Lo.stream = stream;
% % %     end
% % %     
% % % fclose(outfile);
% % % Lo.file = sub_fin;
% % % Lo
% 
% stream = sub_fin(1:1+7,:,:);
% for i=2:height/8   
%     tmp = sub_fin((i-1)*8+1:i*8,:,:);    
%     stream = cat(2,stream,tmp);
% end 
% L1 = (Beginrow/8) * width + Begincol + 1;
% tstream = stream(:, L1: L1+Length-1, :);
% Lo.stream = tstream;