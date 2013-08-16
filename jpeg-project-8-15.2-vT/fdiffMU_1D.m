function diff= fdiffMU_1D(midR, upR)
% endUp = double (upR(end,:)');
% beginMid = double ( midR(1,:)');

endUp = double (upR);
beginMid = double ( midR);

% % % diff = sum(sum(abs(endUp - beginMid)));

powerM=(endUp-beginMid).^2;
sumM = powerM(:,:,1) + powerM(:,:,2) + powerM(:,:,3);
sqrtM = sqrt(sumM);
diff=sum(sqrtM);
