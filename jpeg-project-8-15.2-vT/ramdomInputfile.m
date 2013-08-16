function ramdomInputfile(Infilename,Nocovered,BLOCK, numbers,Rfilename)
% clc;clear;fclose('all');
% Infilename = 'pic_0174.jpg';
% Nocovered = 5;
% BLOCK = 4096;
% numbers = 6;
% Rfilename = 'r_0174.jpg';

orders = frandom_seq(numbers)
% orders = [1,2,3,4,5];
FIDIn  = fopen( Infilename,'r');
fseek(FIDIn,0,'eof');
lsize  = ftell(FIDIn);
frewind(FIDIn);
BufferFIDIn = fread(FIDIn, lsize,'uint8=>uint8');
frewind(FIDIn);
fclose(FIDIn);

streamheader = BufferFIDIn(1:(Nocovered*BLOCK));
streamtailer = BufferFIDIn(BLOCK*floor(lsize/BLOCK)+1:end);
sizetailer   = size(streamtailer,1);
sizeleft     = lsize - size(streamheader,1) - sizetailer;
noblocks     = floor(sizeleft/ (numbers * BLOCK));
sizeblock    = noblocks * BLOCK;
sizeblock
for i = 1: (numbers-1)    
    stream(i).stream = BufferFIDIn((Nocovered*BLOCK)+(i-1)*sizeblock +1 : (Nocovered*BLOCK)+i*sizeblock);
end;
i = i+1;
stream(i).stream = BufferFIDIn((Nocovered*BLOCK)+(i-1)*sizeblock +1 : end-sizetailer);

streambuf = streamheader;
for i=1:numbers
%     disp(orders(i));
    streambuf = cat(1, streambuf,stream(orders(i)).stream);
end
% % streambuf = cat(1, streambuf, streamtailer);
% lsize
% size(streambuf)
writeBuf2File(Rfilename,streambuf);
% x= find(streambuf ~= BufferFIDIn);
% size(x)