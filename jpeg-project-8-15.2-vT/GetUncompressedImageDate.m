function [UncompressedImageDate,IndexOfFragmentPoint1] = GetUncompressedImageDate(FileInfo,IndexOfFragmentPoint,Buffer,NumOfCoveredBlock,TotalNumOfBlock)
%%% maintained by xubenben(HIT)

global Beginrow;
global Begincol;
TempImage='Temp.jpg';
CoveredBuf=Buffer(1).stream;
height=FileInfo.Height;
wid=FileInfo.Width;
SecurityWid= ceil(  ( TotalNumOfBlock/(height/8) )*1.5  );
SecurityWid1= ceil( TotalNumOfBlock/(height/8) );

% L is the number of fragmentPoints
L=size(IndexOfFragmentPoint,2);

%set the interval
for k=1:1:L
    IndexOfFragmentPoint1(2*k-1)=floor((IndexOfFragmentPoint(k)/height)*TotalNumOfBlock-2);
    IndexOfFragmentPoint1(2*k)=floor((IndexOfFragmentPoint(k)/height)*TotalNumOfBlock+2);
end

%initial
for i=2:IndexOfFragmentPoint1(1)
    tmp = Buffer(i).stream;
    CoveredBuf = cat(1, CoveredBuf, tmp);    %cat(1,..), cat matrix function
end;

writeBuf2File(TempImage,CoveredBuf);
EndLocation = FindLocationRecovered(TempImage);
Beginrow=EndLocation.row;
Begincol=EndLocation.col;


%% adjust to confirm that the fragmentpoint is already inthe interval
% IndexOfFragmentPoint1(2*L-1)=IndexOfFragmentPoint1(2*L-1)-200;
% IndexOfFragmentPoint1(2*L)=IndexOfFragmentPoint1(2*L)-200;
%adjust the floor of intervals
k=1;

EndLocation1=0;
EndLocation=0;

while k<=L
    
    CoveredBuf=Buffer(1).stream;
    
    for i=2:IndexOfFragmentPoint1(2*k-1)
        tmp = Buffer(i).stream;
        CoveredBuf = cat(1, CoveredBuf, tmp);    %cat(1,..), cat matrix function
    end;
    EndLocation1=EndLocation;
    
    writeBuf2File(TempImage,CoveredBuf);
    EndLocation = FindLocation(TempImage);
    
    
    if(EndLocation.row>=IndexOfFragmentPoint(k))
        
        tmp= max( floor( ( ( EndLocation.row-IndexOfFragmentPoint(k) )/8  )*SecurityWid1 )/2,SecurityWid1);
        IndexOfFragmentPoint1(2*k-1)=IndexOfFragmentPoint1(2*k-1)-tmp;
        
        IndexOfFragmentPoint1(2*k)=IndexOfFragmentPoint1(2*k-1)+SecurityWid1;
        
        if k==L && EndLocation.row == EndLocation1.row && EndLocation.col == EndLocation1.col
            IndexOfFragmentPoint1(2*k-1)=floor( (IndexOfFragmentPoint1(2*k-1)+IndexOfFragmentPoint1(2*k-3))/2 );
            IndexOfFragmentPoint1(2*k)=IndexOfFragmentPoint1(2*k-1)+SecurityWid1;
        end
        
    else
        k=k+1;
        Beginrow=EndLocation.row;
        Begincol=EndLocation.col;
    end
end


%adjust the edge
if(NumOfCoveredBlock>IndexOfFragmentPoint1(1))
    IndexOfFragmentPoint1(1)=NumOfCoveredBlock;
end

if(TotalNumOfBlock<IndexOfFragmentPoint1(2*L))
    IndexOfFragmentPoint1(2*L)=TotalNumOfBlock;
end

%adjust the ceil of intervals
k=1;
while k<=L
    
    CoveredBuf=Buffer(1).stream;
    
    for i=2:IndexOfFragmentPoint1(2*k)
        tmp = Buffer(i).stream;
        CoveredBuf = cat(1, CoveredBuf, tmp);    %cat(1,..), cat matrix function
    end;
    
    EndLocation1=EndLocation;
    
    
    writeBuf2File(TempImage,CoveredBuf);
    EndLocation = FindLocation(TempImage);
    
    if(EndLocation.row<=IndexOfFragmentPoint(k))
        
        tmp= max( floor(  ( ( ( IndexOfFragmentPoint(k) -EndLocation.row )/8  ) * SecurityWid1 )/2  ),SecurityWid1);
        
        IndexOfFragmentPoint1(2*k)=IndexOfFragmentPoint1(2*k)+tmp;
        
        IndexOfFragmentPoint1(2*k-1)=IndexOfFragmentPoint1(2*k)-SecurityWid1;
        
        if EndLocation.row==IndexOfFragmentPoint(k)+1
            IndexOfFragmentPoint1(2*k)=IndexOfFragmentPoint1(2*k)+SecurityWid1;
        end
        
        if k==L && EndLocation.row == EndLocation1.row && EndLocation.col == EndLocation1.col
            IndexOfFragmentPoint1(2*k-1)=floor( (IndexOfFragmentPoint1(2*k-1)+IndexOfFragmentPoint1(2*k-3))/2 );
            IndexOfFragmentPoint1(2*k)=IndexOfFragmentPoint1(2*k-1)+SecurityWid1;
        end
        
        IndexOfFragmentPoint1(2*L)= min( IndexOfFragmentPoint1(2*L),TotalNumOfBlock);
        
    else
        k=k+1;
        Beginrow=EndLocation.row;
        Begincol=EndLocation.col;
    end
end

tempp=0;
j=1;
flag=1;
% deal with trad
for k=1:L
    
    if flag==1
        tempp(2*j-1)=IndexOfFragmentPoint1(2*k-1);
    end
    
    if k~=L && IndexOfFragmentPoint1(2*k)>=IndexOfFragmentPoint1(2*k+1)
        flag=0;
        continue;
        
    else
        tempp(2*j)=IndexOfFragmentPoint1(2*k);
        j=j+1;
        flag=1;
    end
end
L=j-1;
IndexOfFragmentPoint1=tempp;



for k=1:L
    IndexOfFragmentPoint1(2*k)=IndexOfFragmentPoint1(2*k)+SecurityWid;
end

if(TotalNumOfBlock<IndexOfFragmentPoint1(2*L))
    IndexOfFragmentPoint1(2*L)=TotalNumOfBlock;
end

%% Uncompress the bocks in the intervals!


CoveredBuf=Buffer(1).stream;
i=0;

for i=2:IndexOfFragmentPoint1(1)
    tmp = Buffer(i).stream;
    CoveredBuf = cat(1, CoveredBuf, tmp);    %cat(1,..), cat matrix function
end;

writeBuf2File(TempImage,CoveredBuf);
EndLocation = FindLocationRecovered(TempImage);
UncompressedImageDate(i)=EndLocation;
Beginrow=EndLocation.row;
Begincol=EndLocation.col;



for k=1:1:L
    for i=IndexOfFragmentPoint1(2*k-1)+1:IndexOfFragmentPoint1(2*k)-1
        
        tmp = Buffer(i).stream;
        CoveredBuf = cat(1, CoveredBuf, tmp);
        writeBuf2File(TempImage,CoveredBuf);
        
        EndLocation = FindLocation(TempImage);
        
        Beginrow=EndLocation.row;
        Begincol=EndLocation.col;
        
        UncompressedImageDate(i)=EndLocation;
        
        disp(EndLocation);
        
    end;
    
    
    if k~=L
        for p=IndexOfFragmentPoint1(2*k):IndexOfFragmentPoint1(2*k+1)
            tmp = Buffer(p).stream;
            CoveredBuf = cat(1, CoveredBuf, tmp);
        end
    else
        for p=IndexOfFragmentPoint1(2*k):TotalNumOfBlock
            tmp = Buffer(p).stream;
            CoveredBuf = cat(1, CoveredBuf, tmp);
        end
    end
end
end

