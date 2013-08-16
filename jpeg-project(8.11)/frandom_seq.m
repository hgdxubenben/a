function zorder = frandom_seq(numbers)
a=numbers;
zorder=[1,3,4,2,5,6];
% a = random ( 'uniform', 1, numbers, 5 );
% a = floor(a);
% zorder = zeros(1,numbers)+numbers;
% t = 1;
% flag = 1;
% for i = 1 : 5 
%     for j = 1 : 5  
%         for z = 1 : numbers
%             if a(i,j) == zorder( z )
%                 flag = 1;
%                 break;
%             else
%                 flag = 0;
%             end
%         end
%         if flag == 0
%             zorder(t) = a(i,j);
%             t = t+1;
%         end
%     end;
% end;
% Mod1 = zeros(1,numbers);
% for i = 1 : numbers    
%     if zorder(i) ~= numbers
%         Mod1(zorder(i)) = 1;
%     end    
% end
% 
% t = numbers;
% for i = 1 : (numbers - 1)
%     if Mod1(i) == 0
%         zorder(t) = i;
%         t = t-1;
%     end
% end;