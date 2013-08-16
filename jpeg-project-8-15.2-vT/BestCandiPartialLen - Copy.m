function Ibc = BestCandiPartialLen_copy(upR, nblocks, stm, width, Znum)

value = zeros(Znum,1) + 99999999;
vector = zeros(Znum,1) - 1;
[Maxvalue, Locvector] = max(value);

tmptc(1:nblocks)=zeros + 99999999;
I1 = -1; I2 = -1;
minI1 = 99999999;
minI2 = 99999999;

%upR = upR(:,end-width+1:end,:);
% * * * * Mid stream
for x = 1 : nblocks
    if stm(x).order == -1 && stm(x).order ~= -2
%         disp('....');
%         Ttm = Ttm + 1;
        
        szX = size(stm(1,x).stream,2); 
        stmX = stm(x).stream;
        if szX > width
            stmX = stm(x).stream(:, 1:width, :);
            szX = width;
        end
        
%%        
        szT = floor(6*szX/9);
        if szT <= width            
            midR = stmX(1,1:szT,:);
            upRt = upR(end,1:szT,:);
            MUT1 = fdiffMU_1D(midR,upRt);
            if MUT1 > Maxvalue * szX 
%                  Ttm1 = Ttm1 + 1;
                continue;
            end
        else            
            MU1 = 0;
        end;
        
        szT2 = floor(7*szX/9);  
        if szT2 <= width
            midR = stmX(1,szT+1:szT2,:);
            upRt = upR(end,szT+1:szT2,:);
            MUT2 = MUT1+fdiffMU_1D(midR,upRt);
            if MUT2 > Maxvalue * szX 
%                  Ttm2 = Ttm2 + 1;
                continue;
            end
        else
            MUT2 = MU1;
        end
        
        szT3 = floor(8*szX/9);           
        if szT3 <= width
            midR = stmX(1,szT2+1:szT3,:);
            upRt = upR(end,szT2+1:szT3,:);
            MUT3 = MUT2+fdiffMU_1D(midR,upRt);
            if MUT3 > Maxvalue * szX 
%                  Ttm3 = Ttm3 + 1;
                continue;
            end
        else
            MUT3 = MUT2;
        end;
        
%%                
        midR = stmX(1,szT3+1:end,:);
        upRt = upR(end,szT3+1:szX,:);		
        
		MUT4 = fdiffMU_1D(midR,upRt);
		TmptC = (MUT3+MUT4) / szX;        

        tmptc(x) = TmptC;        
%         Ttm4 = Ttm4 + 1;
%         disp(TmptC);
%         disp(size(stmX)); 
%         disp(size(upR));
        %%
        if TmptC < Maxvalue
            value(Locvector) = TmptC;
            vector(Locvector) = x;
            [Maxvalue, Locvector] = max(value);
%             [Maxvalue, Locvector] = max(value);
        end        
        
   end
end
% % % disp('vector');
% % % disp(vector);
[minI1, m] = min(value);
value(m) = -10;
[minI2, t] = max(value);
I1 = vector(m);
I2 = vector(t);
if I2 == -1
    I2 = I1;
    minI2 = minI1;
end
Ibc(1) = I1;
Ibc(2) = I2;
Ibc(3) = minI1;
Ibc(4) = minI2;
% disp(tmptc');
% Ibc.value = value;
% Ibc
% [minI1, m] = min(value);
% [minI3, t] = max(value);
% I1 = vector(m);
% I3 = vector(t);
% if I3 == -1
%     I3 = I1;
%     minI3 = minI1;
% end
% value(m) = 99999999;
% [minI2, s] = min(value);
% I2 = vector(s);
% if I2 == -1
%     I2 = I1;
%     minI2 = minI1;
% end
% 
% Ibc(1) = I1;
% Ibc(2) = I2;
% Ibc(3) = I3;
% Ibc(4) = minI1;
% Ibc(5) = minI2;
% Ibc(6) = minI3;