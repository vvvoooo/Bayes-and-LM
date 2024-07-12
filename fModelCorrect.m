function [Y] = fModelCorrect(p)
fid=fopen('bianliang.txt','wt');  %导入初始变量参数
fprintf(fid,'%16.8f%16.8f%16.8f',p);
fclose(fid);
system(strcat('"F:\ANSYS Inc\v221\ANSYS\bin\winx64\ANSYS221.exe"',  ' -b -p ane3fl -i', ' F:\guopeiyu\zhuiqiao\012\banapdlcorrect.txt', ' -o', ' F:\guopeiyu\zhuiqiao\012\process.out'));
Y = load('celiang1.txt');          
%T2 = load('result2.txt'); % measurement vector
%T1 = [T1 
%   ];
end

