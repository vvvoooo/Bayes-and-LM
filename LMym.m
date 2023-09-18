% Levenberg-Marquardt Method 三维结构

p = [1;1;1]; % Guess parameters
fid=fopen('bianliang.txt','wt');  %导入初始变量参数
fprintf(fid,'%8.4f%8.4f%8.4f',p);
fclose(fid);
N = size(p, 1);
I = 19;
%打开目标ansys并带入参数进行运算
system(strcat('"F:\ANSYS Inc\v221\ANSYS\bin\winx64\ANSYS221.exe"',  ' -b -p ane3fl -i', ' F:\guopeiyu\apdl\5\banapdl.txt', ' -o', ' F:\guopeiyu\apdl\5\process.out'));
Y = load('celiang.txt');          % measurement vector
delp = 10^(-2);
J = zeros(I, N);            % sensitivity matrix 
e1 = 3;
e2 = 1.5;
e3 = 3;
mu = 0.0001;
% plot(1:200, Y);

count = 0;
while true
 count = count + 1; 
 %ansys得到初始结果
 system(strcat('"F:\ANSYS Inc\v221\ANSYS\bin\winx64\ANSYS221.exe"',  ' -b -p ane3fl -i', ' F:\guopeiyu\apdl\5\banapdl.txt', ' -o', ' F:\guopeiyu\apdl\5\process.out'));
 T = load('result.txt');    %读取初始温度矩阵
 S = sum((Y-T).^2);         % Squared Error 
 
 Tpj = load('result.txt');                     %读取未改变的参数带来的温度矩阵
 
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         pInc(j) = pInc(j) + delp*pInc(j);
         
         fid=fopen('bianliang.txt','wt');              %每次循环中的参数导入过程
         fprintf(fid,'%8.4f%8.4f%8.4f',pInc);
         fclose(fid);
         %打开ansys运行得到新的温度矩阵
         system(strcat('"F:\ANSYS Inc\v221\ANSYS\bin\winx64\ANSYS221.exe"',  ' -b -p ane3fl -i', ' F:\guopeiyu\apdl\5\banapdl.txt', ' -o', ' F:\guopeiyu\apdl\5\process.out'));
         Tdelpj = load('result.txt');    
  
         J(i,j) = (Tdelpj(i) - Tpj(i))/(delp*p(j));
     end
 end

 omega = diag(diag(J'*J));   

 pupdt = p;

 while true
     
     fModelGuess = load('result.txt');
     pupdt = pupdt + (pinv(J'*J + mu*omega))*(J')*(Y - fModelGuess);   
     % updated T
     fid=fopen('bianliang.txt','wt');  %导入变化变量参数
     fprintf(fid,'%8.4f%8.4f%8.4f',pupdt);
     fclose(fid);
     %打开ansys运行得到新的温度矩阵
     system(strcat('"F:\ANSYS Inc\v221\ANSYS\bin\winx64\ANSYS221.exe"',  ' -b -p ane3fl -i', ' F:\guopeiyu\apdl\5\banapdl.txt', ' -o', ' F:\guopeiyu\apdl\5\process.out'));
     T = load('result.txt');


     Supdt = sum((Y-T).^2);
     if Supdt >= S         
       mu = 10*mu;
     else
       mu = 0.1*mu;  
       break;
     end
 end

 if (Supdt<e1 || norm(J'*(Y - fModelGuess))<e2 || norm(pupdt-p)<e3)
     pEstimate = pupdt;
     break;
 else
     p = pupdt;
 end
 
end

% t = linspace(0,40,200);
% qt=30.*sin(t);
% qtyuce=cos(t).*pEstimate(1)+sin(t).*pEstimate(2)+ t.^2.*pEstimate(3);
% plot(t,qt,'-r',t,qtyuce,'-xb');
