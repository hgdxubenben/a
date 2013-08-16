function writeBuf2File(Outfilename, buff)
global writetime;
btime = cputime;
outfile = fopen(Outfilename,'wb'); 
while outfile == -1
    outfile = fopen(Outfilename,'wb');        
end
fwrite(outfile,buff);
% fclose(outfile);
fclose('all');
etime = cputime - btime;
writetime = writetime + etime;