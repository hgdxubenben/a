function diff=fdiffML_1D(midR,leftR)
% endL=double(leftR(:,end,:));
% beginM = double(midR(:,1,:));
endL = double(leftR);
beginM = double(midR);

% % % diff = sum(sum(abs(endUp - beginMid)));

powerM=(endL-beginM).^2;
sumM = powerM(:,:,1) + powerM(:,:,2) + powerM(:,:,3);
sqrtM = sqrt (sumM);
diff=sum(sqrtM);