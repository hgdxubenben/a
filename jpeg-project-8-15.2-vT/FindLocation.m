function Lo=FindLocation(filename)

global Beginrow;
global Begincol;
global timefind;
global Lotimes;
global sample_factor_v;
global sample_factor_h;
btime = cputime;
Lotimes = Lotimes + 1;
% -- find the end pixel location

outfile = fopen(filename);  
while outfile == -1
    outfile = fopen(filename);        
end



sub_fin = imread(filename);
width = size(sub_fin, 2);
height = size(sub_fin, 1);

if (mod(width,sample_factor_v*8) ~= 0)
    width = width - mod(width, sample_factor_v*8);
    sub_fin = sub_fin(:,1:width,:);
end;
if (mod(height,sample_factor_h*8) ~= 0)
    height = height - mod(height, sample_factor_h*8);
    sub_fin = sub_fin(1:height,:,:);
end;

% sample_factor_v = 2;
% sample_factor_h = 2;
% sub = sub_fin(:,:,1);

row_no = flipud(sub_fin(:,end,1));
row    = find(row_no ~= 128, 1);
row    = height - row + 1;
if mod(row,8*sample_factor_h)~= 0 
    row = round(row/(sample_factor_h*8)) * (sample_factor_h*8);
end;

if row == height 
    column = width;
else
    col_no = fliplr(sub_fin(row+3,:,1));
    column = find(col_no ~= 128, 1);
    if sum(column) == 0 
        disp('sum column = 0');    
        column = width;
    end
    column = width - column + 1;
    column = mod((column - sample_factor_v*8 + width),width);
end
column = column - sample_factor_v*8;
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
row = row + 1;
column = column + 1;

Length = floor((row - Beginrow)/(sample_factor_h*8) - 1)*width + (width-Begincol-1) + column;
Lo.len = Length;



% ---- combine and find the data stream
% * * * * combine
    if Lo.row ~= Beginrow
        stream = sub_fin( Beginrow+1 : Beginrow+sample_factor_h*8, Begincol+1 : width,:);
        for i= (Beginrow/(8*sample_factor_h) +2): (row/(sample_factor_h*8))
            tmp = sub_fin((i-1)*(sample_factor_h*8)+1 : i*(sample_factor_h*8),:,:);
            stream = cat(2,stream, tmp);
        end
        if Lo.row ~= height
            tmp = sub_fin(Lo.row+1: Lo.row+(sample_factor_h*8), 1:Lo.col, :);
            stream = cat(2, stream, tmp);
        else
            tmp = sub_fin(Lo.row-(sample_factor_h*8)+1 : Lo.row, 1:Lo.col,:);
            stream = cat(2, stream, tmp);
        end
        Lo.stream = stream;
    else
        stream = sub_fin( Beginrow+1 : Beginrow + sample_factor_h*8, Begincol+1 : Lo.col,:);
        Lo.stream = stream;
    end
    
fclose(outfile);
etime = cputime - btime;
timefind = timefind + etime;