function [T] = fModelGuess(p)
global ps1 ps2 ps3 ps4 NN;
if NN==100;
 fid=fopen('bianliang.txt','wt');  %导入初始变量参数
 fprintf(fid,'%16.8f%16.8f%16.8f',p);
 fclose(fid);
else 
 fid=fopen('bianliang.txt','wt');  %导入初始变量参数
 fprintf(fid,'%16.8f%16.8f%16.8f',ps1);
 fclose(fid);
end
if NN==200;
    fid=fopen('bianliang.txt','wt');  %导入初始变量参数
    fprintf(fid,'%16.8f%16.8f%16.8f',p);
    fclose(fid);
else 
    fid=fopen('bianliang.txt','wt');  %导入初始变量参数
    fprintf(fid,'%16.8f%16.8f%16.8f',ps2);
    fclose(fid);
end
if NN==300;
    fid=fopen('bianliang.txt','wt');  %导入初始变量参数
    fprintf(fid,'%16.8f%16.8f%16.8f',p);
    fclose(fid);
else 
    fid=fopen('bianliang.txt','wt');  %导入初始变量参数
    fprintf(fid,'%16.8f%16.8f%16.8f',ps3);
    fclose(fid);
end
if NN==400;
    fid=fopen('bianliang.txt','wt');  %导入初始变量参数
    fprintf(fid,'%16.8f%16.8f%16.8f',p);
    fclose(fid);
else 
    fid=fopen('bianliang.txt','wt');  %导入初始变量参数
    fprintf(fid,'%16.8f%16.8f%16.8f',ps4);
    fclose(fid);
end
system(strcat('"F:\ANSYS Inc\v221\ANSYS\bin\winx64\ANSYS221.exe"',  ' -b -p ane3fl -i', ' F:\guopeiyu\zhuiqiao\012\banapdl.txt', ' -o', ' F:\guopeiyu\zhuiqiao\012\process.out'));
T = load('result1.txt');          
%T2 = load('result2.txt'); % measurement vector
%T = [T1 
%    ];
end




