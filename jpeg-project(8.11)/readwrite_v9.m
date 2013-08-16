% function readwrite_v9(OriginalFiContactLename, InFiContactLename,OutFiContactLename)
% consider the next BlockSize as an available variable.... updated from v4
% consider the ContactLength of min(Isz, Ibc1, Ibc2)
% consider the identified next BlockSize could be used in next round.
%%  % Initial the program
fclose('all');
clc
clear
warning off;
global Beginrow;    %recovered picture's row & col
global Begincol;
global writetime;    %testing: calculate time comsumption
global timefind;
global Lotimes;
global sample_factor_v;  % MCU'width
global sample_factor_h;  % MCU'hight
sample_factor_v = 2;   %v=2, mcu's v pixel = 2 * 8=16
sample_factor_h = 1;
writetime = 0;
timefind  = 0;
tmsnext   = 0;
Lotimes = 0;
OriginalFiContactLename  = 'DSC07268.jpg';       %input original file name
InFiContactLename   = 'Randomedfile.jpg';
% InFiContactLename   = 'carv15.jpg';
% InFiContactLename   = OriginalFiContactLename;
OutFiContactLename  = 'covered1.jpg';
fiContactLename     = OutFiContactLename;

BlockSize		= 4096;             %data BlockSize size
numbers     = 6;                %number of fragments of Original file
NumOfCoveredBlock	= 15;               %header's blcok no.(init), recovered's BlockSize total numbers.
Init_NumOfCovered = NumOfCoveredBlock;

%random the file fragments
ramdomInputfile(OriginalFiContactLename,NumOfCoveredBlock,BlockSize,numbers,InFiContactLename);

info   = imfinfo(InFiContactLename);
width  = info.Width;
height = info.Height;

FID_IN  = fopen( InFiContactLename,'r');
fseek(FID_IN,0,'eof');
InFileSize  = ftell(FID_IN);
frewind(FID_IN);
BufferFID_IN = fread(FID_IN, InFileSize,'uint8=>uint8');
fclose(FID_IN);

%%  % Organization of the input BlockSize file.
for i = 1: NumOfCoveredBlock                        %preparation of file header's data BlockSize
    covered(i) = i;
    available(i) = -1;
    Buffer(i).available = -1;
    Buffer(i).stream = BufferFID_IN( (i-1)*BlockSize+1 : i*BlockSize);
end
for i = NumOfCoveredBlock+1 : floor(InFileSize/BlockSize)    %preparation of all unrecovered data BlockSize
    available( i) = i;
    Buffer(i).available = i;
    Buffer(i).stream = BufferFID_IN( (i-1)*BlockSize+1 : i*BlockSize);
end
if floor(InFileSize/BlockSize) ~= ceil(InFileSize/BlockSize)  %preparation of last fragment's data BlockSize
    i=i+1;
    Buffer(i).stream = BufferFID_IN((i-1)*BlockSize+1 :end);
    Buffer(i).available = -1;
    available(i) = -1;   %exclude the last block
end;

BufferFID_IN = 0; % Release the Buffer of Input file's data stream

NumOfBlock = i-1;  % sum of all BlockSizes' exclude the last BlockSize

%get the 15 head block
CoveredBuf = Buffer(1).stream;
for i=2:NumOfCoveredBlock
    tmp = Buffer(i).stream;
    CoveredBuf = cat(1, CoveredBuf, tmp);    %cat(1,..), cat matrix function
end;

%%  main function initial ..
Btime = cputime;    %Begin time
Beginrow = 0;
Begincol = 0;
tms = 1;
index = NumOfCoveredBlock;      %recovered BlockSize's no(index)

%the average similarity between initial lines
meanI = MeanInit(fiContactLename, CoveredBuf);

value(1:NumOfCoveredBlock)      = zeros;            % value: store the recovered similarity value of each fragment
value(Init_NumOfCovered)   = meanI;
Inext    = -1;
Ibc1next = -1;
Ibc2next = -1;
tnext    = 0;
TotalNumOfBlock = floor(InFileSize/BlockSize);
%%  % main function
% for BlockIndex=NumOfCoveredBlock+1:60
for BlockIndex=NumOfCoveredBlock+1:TotalNumOfBlock
    
    if (NumOfCoveredBlock - 10) >= Init_NumOfCovered
        meanI = mean(value(end-10:end-1)); %avrage
        fprintf('------NumOfCoveredBlock-10\t%d\n',meanI);
    else
        meanI = mean(value(Init_NumOfCovered:(end-1)));
        fprintf('------NumOfCoveredBlock-10-else\t%d\n',meanI);
    end
    
    fprintf('\n------------------meanI-------BlockIndex-----------\n');
    fprintf('meanI = %2.3d \t\t BlockIndeximes: %d\t\tindex=%d\n', meanI, BlockIndex,index);
    notfornext = 1;
    flagup = 1;
    writeBuf2File(OutFiContactLename, CoveredBuf);
    EndLocation = FindLocationRecovered(OutFiContactLename)    %stream(1).可以使用原?以?算出的信息-improvment point
    stream(1).stream	= EndLocation.stream;
    stream(1).row   	= EndLocation.row;
    stream(1).col		= EndLocation.col;
    stream(1).ContactLen   	= EndLocation.len;
    stream(1).order 	= 1;
    Beginrow 			= EndLocation.row;
    Begincol 			= EndLocation.col;
    
    upR = stream(1).stream;
    tms = 1;
    if 16 == EndLocation.len
        disp('ContactLength == 16'); %  disp(X)  :   If X is a string, the text is displayed.
        break;
    end;
    
    
    switch (index+1)
        case Inext
            value(NumOfCoveredBlock+1) = Imu2;
        case Ibc1next
            value(NumOfCoveredBlock+1) = Imuin1;
        case Ibc2next
            value(NumOfCoveredBlock+1) = Imuin2;
        otherwise
            value(NumOfCoveredBlock+1) = 99999;
    end
    %% index+1 ContactLen=32'stream compare with upR stream -- checking the physically nearby data BlockSize
    % 	if (NumOfCoveredBlock ~= Init_NumOfCovered) && ((index+1) <= TotalNumOfBlock)
    if ((index+1) <= TotalNumOfBlock && Buffer(index+1).available ~= -1)
        disp('after');
        NextBuf 	= cat(1, CoveredBuf, Buffer(index+1).stream);
        writeBuf2File(OutFiContactLename, NextBuf);
        NextLocation 	= FindLocation(OutFiContactLename)
        NextStream = NextLocation.stream;
        
        ContactLen = min(size(NextStream,2),width);
        edi        = fdiffMU_1D(upR(end-1,1:ContactLen,:),    upR(end,1:ContactLen,:))/ContactLen; %edi:upR's last two rows's similarity
        ednear    = fdiffMU_1D(NextStream(1,1:ContactLen,:),upR(end,1:ContactLen,:))/ContactLen; %ednear:similarity between upR's last row and  NextStream's first row
        
        EDi(BlockIndex)  = ednear;
        EDNear(BlockIndex)  = edi;
        EDLen(BlockIndex)= ContactLen;
        ED(BlockIndex)= (ednear - edi);
        
        if ednear <= edi || (ednear - edi) <= (ContactLen*0.007)
            tmsnext = tmsnext +1;
            index                   = index+1;
            available(index)        = -1;
            covered(NumOfCoveredBlock+1)    = index;  % stand for the rcovered block order
            CoveredBuf              = cat(1, CoveredBuf, Buffer(index).stream); % stand for the covered image buffer
            Buffer(index).available = -1;
            NumOfCoveredBlock               = NumOfCoveredBlock +1;
            Inext                   = -1;
            Ibc1next                = -1;
            Ibc2next                = -1;
            ContactLengthmatrix(BlockIndex)        = size(NextStream,2);
            value(BlockIndex)               = ednear;
            continue;
        end;
    end;
    
    fprintf('value(NumOfCoveredBlock+1)=%d\t Buffer(index+1).able=%d\n', value(NumOfCoveredBlock+1),Buffer(index+1).available );
    if value(NumOfCoveredBlock+1) < (meanI + 4) && Buffer(index+1).available ~= -1          % tolerate for dramatically change
        fprintf('\t\t*******Yes for next....%d\t%d',tnext,BlockIndex);
        tnext = tnext+1;
        index                   = index+1;
        available(index)        = -1;
        covered(NumOfCoveredBlock+1)    = index;
        CoveredBuf              = cat(1, CoveredBuf, Buffer(index).stream);
        Buffer(index).available = -1;
        NumOfCoveredBlock               = NumOfCoveredBlock +1;
        Inext                   = -1;
        Ibc1next                = -1;
        Ibc2next                = -1;
        fornextflag = 1;
        continue;
    else
        fornextflag = -1;
    end
    
    %% fragmentpoint
    Couldbenext = -1;
    Incnext     = -1;
    if CoveredBuf(end) == 255 && (Buffer(index+1).stream(1) ~= 0)  %255=0xFF, JPEG data stream: 0xFF - 0x00. otherwise D8,D9, 0xFFD8
        CouldbeDone = -1;
        disp('does it not functional!!!');
    else
        CouldbeDone = 1;
        disp('It Does functional!!!');
    end;
    %% -- -- -- next BlockSize..
    %I delete this code as they will not ues
%     if CouldbeDone == 1 && (index+1) <= NumOfBlock && (Buffer(index+1).available ~= -1)
%         Couldbenext = 1;
%         nextBlock = Buffer(index+1).stream;
%         in3 = cat(1, CoveredBuf, nextBlock);
%         writeBuf2File(OutFiContactLename, in3);
%         nextLocation = FindLocation(OutFiContactLename);
%         nextStream = nextLocation.stream;
%         Isz = min(size(nextStream,2), width);
%         Imu = fdiffMU_1D(nextStream(1,1:Isz,:),upR(end,1:Isz,:))/Isz;
%         Istream = nextStream;
%         IszContactLen = Isz;
%         fprintf('------Isz & Imu------: Isz = %d\t\tImu=%d\n',Isz, Imu);
%         
%         if (Imu <=  meanI*3 ) &&((index + 2) <= NumOfBlock) && (Isz < (2/3 * width)) && Buffer(index+2).available ~= -1
%             if Buffer(index+1).stream(end) == 255 && Buffer(index+2).stream(1) ~= 0
%                 disp('Next next stream 255 !0');
%             else
%                 Inext = index + 2;
%                 nextBlockSize2 = Buffer(index+2).stream;
%                 in3 = cat(1, in3, nextBlockSize2);
%                 writeBuf2File(OutFiContactLename, in3);
%                 nextLocation2 = FindLocation(OutFiContactLename);
%                 nextStream2 = nextLocation2.stream;
%                 Isz2 = min(size(nextStream2,2), width);
%                 Imu2 = fdiffMU_1D(nextStream2(1,Isz+1:Isz2,:),upR(end,Isz+1:Isz2,:))/(Isz2-Isz);
%                 %             Imu2 = (fdiffMU_1D(nextStream2(1,Isz+1:Isz2,:),upR(end,Isz+1:Isz2,:)) + Imu*Isz)/Isz2;
%                 %             fprintf('------Isz & Imu-next-: Isz2=%d \t\t \n', Isz2);
%                 Imu = (Imu*Isz+Imu2*(Isz2-Isz))/Isz2;
%                 fprintf('------Isz & Imu-next-: Isz2=%d \t\t Imu2=%d\BlockIndexvalue=%d\n', Isz2, Imu2,Imu);
%                 Istream = nextStream2;
%                 IszContactLen = Isz2;
%                 in3 = 0;
%             end
%         end
%         fprintf('* * Imu=%d,,, meanI=%d\n', Imu, meanI);
%         if  (Imu <= meanI+4) %&& (-1 ~= fornextflag)
%             fprintf(' * * Imu=%d,,, meanI=%d\n', Imu, meanI);
%             Incnext = 1;
%             index = index + 1;
%             disp('Imu <= meanI+4** ');
%         end
%     end
    %% %searching for all avialable candidate BlockSizes
    if Incnext == -1   
        for i=1 : NumOfBlock
            if Buffer(i).available == -1
                continue;
            end;
            if CoveredBuf(end) == 255 && Buffer(available(i)).stream(1) ~= 0
                disp('does it not functional?>>');
                continue;
            end;
            
            in2 = Buffer(available(i)).stream;
            in3 = cat(1, CoveredBuf, in2);
            writeBuf2File(OutFiContactLename, in3);
            EndLocation = FindLocation(OutFiContactLename);
            if EndLocation.len <= 16
                continue;
            end;
            tms = tms + 1;
            stream(tms).stream	= EndLocation.stream;
            stream(tms).row		= EndLocation.row;
            stream(tms).col   	= EndLocation.col;
            stream(tms).ContactLen		= EndLocation.len;
            stream(tms).order	= -1;
            stream(tms).Index	= i;
            
            in2 = 0;in3 = 0;
        end
        
        if tms == 1
            fprintf('****tms == 1****\n');
            break;
        end;
        
        Ibc = BestCandiPartialLen(upR, tms, stream, width, 2);   %Ibc: top m candidates
        fprintf('Ibc values: %d\t%d\t%d\t%d\n',Ibc(1),Ibc(2),Ibc(3),Ibc(4));
        %% next BlockSize of Ibc(1) & Ibc(2)     判?BlockSize x, and then merged with top m candidates as II
        Ibc1stream = stream(Ibc(1)).stream;
        Ibc2stream = stream(Ibc(2)).stream;
        Ibc1sz     = size(Ibc1stream,2);
        Ibc2sz     = size(Ibc2stream,2);
        
        index1 = stream(Ibc(1)).Index;
        index2 = stream(Ibc(2)).Index;
        value1 = Ibc(3);
        value2 = Ibc(4);
        fprintf('==Ibc1sz=%d\t\tIbc2sz=%d\n',Ibc1sz, Ibc2sz);
        if stream(Ibc(1)).ContactLen <= 16
            disp('stream(Ibc(1)).ContactLen <= 16');
            break;
        end;
        
        if stream(Ibc(1)).ContactLen <= (2/3 * width) && ((index1 +1) <= NumOfBlock) &&(Buffer(index1 +1).available ~= -1)
            disp('==========damed if ');
            Ibc1next = index1+1;
            in1     = Buffer(index1).stream;
            in1next = Buffer(index1 +1).stream;
            in1buf  = cat(1, in1, in1next);
            in1buf  = cat(1, CoveredBuf, in1buf);
            writeBuf2File(OutFiContactLename, in1buf);
            nextLocation2 = FindLocation(OutFiContactLename);
            in1nextstream = nextLocation2.stream;
            Iszin1 = min(size(in1nextstream,2), width);
            Imuin1 = (fdiffMU_1D(in1nextstream(1,Ibc1sz+1 : Iszin1,:), upR(end,Ibc1sz+1 : Iszin1,:)))/(Iszin1-Ibc1sz);
            value1 = (Imuin1*(Iszin1-Ibc1sz) + value1*Ibc1sz)/Iszin1;
            %             Imuin1 = (fdiffMU_1D(in1nextstream(1,Ibc1sz+1 : Iszin1,:), upR(end,Ibc1sz+1 : Iszin1,:)) + value1*Ibc1sz)/Iszin1;
            fprintf('------Ibc(1)...Iszin1 & Imuin1---%d\t%d-\BlockIndexvalue=%d-\n', Iszin1,Imuin1,value1);
            Ibc1stream = in1nextstream;
            Ibc1sz     = size(in1nextstream,2);
            fprintf('Ibs1sz = %d, sizeofIbc1stream=%d\n',Ibc1sz,size(Ibc1stream,2));
            in1 			= 0;
            in1next 		= 0;
            in1buf  		= 0;
            nextLocation2 	= 0;
            in1nextstream	= [];
        end
        if stream(Ibc(2)).ContactLen <= (2/3 * width) && ((index2 +1) <= NumOfBlock) && (Buffer(index2 +1).available ~= -1)
            Ibc2next = index2+1;
            in1     = Buffer(index2).stream;
            in1next = Buffer(index2 +1).stream;
            in1buf  = cat(1, in1, in1next);
            in1buf  = cat(1, CoveredBuf, in1buf);
            writeBuf2File(OutFiContactLename,in1buf);
            nextLocation2 = FindLocation(OutFiContactLename);
            in2nextstream = nextLocation2.stream;
            Iszin2 = min(size(in2nextstream,2), width);
            Imuin2 = (fdiffMU_1D(in2nextstream(1,Ibc2sz+1:Iszin2,:), upR(end,Ibc2sz+1:Iszin2,:)))/(Iszin2-Ibc2sz);
            value2 = (Imuin2*(Iszin2-Ibc2sz) + value2*Ibc2sz)/Iszin2;
            %             Imuin2 = (fdiffMU_1D(in2nextstream(1,Ibc2sz+1:Iszin2,:), upR(end,Ibc2sz+1:Iszin2,:))+value2*Ibc2sz)/Iszin2;
            fprintf('----Ibc(2)...Iszin2 & Imuin2---%d\t%d\BlockIndexvalue=%d-\n',Iszin2, Imuin2, value2);
            Ibc2stream = in2nextstream;
            Ibc2sz     = size(Ibc2stream,2);
            fprintf('Ibs2sz = %d, sizeofIbc2stream=%d\n',Ibc2sz,size(Ibc2stream,2));
            in1 			= 0;
            in1next 		= 0;
            in1buf  		= 0;
            nextLocation2 	= 0;
            in1nextstream	= [];
        end
        %%  compare the Imu, value1, value2,   %find the best one from II
        flagofindex = -1;
%         if Couldbenext == 1
%             fprintf('if Could be Done == 1\n');
%             if (index+1) == stream(Ibc(1)).Index || ((index+1) == stream(Ibc(2)).Index && abs(value2-value1) <= 20 )
%                 fprintf('00 if index of next happens in Ibc\n');
%                 index = index + 1;
%                 value(NumOfCoveredBlock +1) = Imu;
%                 vector(NumOfCoveredBlock+1) = index;
%                 ContactLengthmatrix(BlockIndex)    = Isz;
%                 flagofindex = 1;
%                 Buffer(index).available = -1;
%             end;
%             if flagofindex == -1
%                 %sz = min(IszContactLen, Ibc1sz);
%                 sz=Ibc1sz;
%                 sz = min(sz,     Ibc2sz);
%                 sz = min(sz,     width);
%                 %  mu1 = fdiffMU_1D(Istream(1,1:sz,:),     upR(end,1:sz,:))/sz;
%                 mu2 = fdiffMU_1D(Ibc1stream(1,1:sz,:),  upR(end,1:sz,:))/sz;
%                 mu3 = fdiffMU_1D(Ibc2stream(1,1:sz,:),  upR(end,1:sz,:))/sz;
%                 if mu2<=mu3
%                     index = stream(Ibc(1)).Index;
%                     value(NumOfCoveredBlock +1) = Ibc(3);
%                     vector(NumOfCoveredBlock+1) = index;
%                     ContactLengthmatrix(BlockIndex)    = stream(Ibc(1)).ContactLen;
%                     Buffer(index).available = -1;
%                 elseif
%                      index = stream(Ibc(2)).Index;
%                     value(NumOfCoveredBlock +1) = Ibc(4);
%                     vector(NumOfCoveredBlock+1) = index;
%                     ContactLengthmatrix(BlockIndex)    = stream(Ibc(2)).ContactLen;
%                     Buffer(index).available = -1;
%                 end
% %                 if mu1 <= mu2 && mu1 <= mu3
% %                     fprintf(' 01 if Imuc is the min...');
% %                     index = index + 1;
% %                     value(NumOfCoveredBlock +1) = Imu;
% %                     vector(NumOfCoveredBlock+1) = index;
% %                     ContactLengthmatrix(BlockIndex)    = Isz;
% %                     Buffer(index).available = -1;
% %                 elseif mu2 <= mu1 && mu2 <= mu3
% %                     disp('if value1 is the min..');
% %                     index = stream(Ibc(1)).Index;
% %                     value(NumOfCoveredBlock +1) = Ibc(3);
% %                     vector(NumOfCoveredBlock+1) = index;
% %                     ContactLengthmatrix(BlockIndex)    = stream(Ibc(1)).ContactLen;
% %                     Buffer(index).available = -1;
% %                 else
% %                     disp('if value2 is the min..');
% %                     index = stream(Ibc(2)).Index;
% %                     value(NumOfCoveredBlock +1) = Ibc(4);
% %                     vector(NumOfCoveredBlock+1) = index;
% %                     ContactLengthmatrix(BlockIndex)    = stream(Ibc(2)).ContactLen;
% %                     Buffer(index).available = -1;
% %                 end
%             end
%         end
        if Couldbenext == -1
            sz = min(Ibc1sz, Ibc2sz);
            sz = min(sz,     width);
            mu2 = fdiffMU_1D(Ibc1stream(1,1:sz,:),  upR(end,1:sz,:))/sz;
            mu3 = fdiffMU_1D(Ibc2stream(1,1:sz,:),  upR(end,1:sz,:))/sz;
            if mu2 <= mu3
                fprintf('if Ibc(1) is the min..');
                index = stream(Ibc(1)).Index;
                value(NumOfCoveredBlock +1) = Ibc(3);
                vector(NumOfCoveredBlock+1) = index;
                ContactLengthmatrix(BlockIndex)    = stream(Ibc(1)).ContactLen;
                Buffer(index).available = -1;
            else
                fprintf('if Ibc(2) is the min..');
                index = stream(Ibc(2)).Index;
                value(NumOfCoveredBlock +1) = Ibc(4);
                vector(NumOfCoveredBlock+1) = index;
                ContactLengthmatrix(BlockIndex)    = stream(Ibc(2)).ContactLen;
                Buffer(index).available = -1;
            end
        end
        
    else
        value (NumOfCoveredBlock+1) = Imu;
        vector(NumOfCoveredBlock+1) = index+1;
        ContactLengthmatrix(BlockIndex) 	= Isz;
    end;
    %     if tms ~= 1
    available(index) 		= -1;
    covered(NumOfCoveredBlock+1)	= index; %Ibc(1)+NumOfCoveredBlock - 1;
    CoveredBuf				= cat(1, CoveredBuf, Buffer(index).stream);
    Buffer(index).available = -1;
    
    fprintf('\t\t\t\tindex =    %d\n',index);
    writeBuf2File(OutFiContactLename, CoveredBuf);
    NumOfCoveredBlock = NumOfCoveredBlock + 1;
    %     end
    %     stream = [];
    %     Ibc
end
break;
Etime = cputime - Btime;
CoveredBuf = cat(1, CoveredBuf, Buffer(end).stream);
writeBuf2File(OutFiContactLename, CoveredBuf);
fprintf('Total Time = %d, writetime=%d\t, timefind = %d\n',Etime, writetime, timefind);
fprintf('Recover Time = %d\n', Etime-writetime-timefind);
fprintf('Running times of find location: %d\n',Lotimes);
fclose(fid);

