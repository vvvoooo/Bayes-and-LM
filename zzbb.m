
%测量数据的长度，或者说单个测点的测量时间个数
p=[300,500,10000];
Y = fModelCorrect(p);
Y100 = Y(1:100);
Y200 = Y(1:200);
Y300 = Y(1:300);
Y400 = Y(1:4000);
I100 = siez(Y100,1);
I200 = siez(Y200,1);
I300 = siez(Y300,1);
I400 = siez(Y400,1);

global ps1 ps2 ps3 ps4 NN;

ps1=[1,1,1];
ps2=[1,1,1];
ps3=[1,1,1];
ps4=[1,1,1];

NN = 100;%该部分代码的计算时间序列为0-100s
Y = Y100;
number = [-100 3 300;-10 5 50;-10 -1 0.1;];%待反演参数的取值范围
m = size(number,1);%待反演参数数量
a = size(number,2);
p = zeros(m,a^m);
%先创造参数空间，在这个参数空间里，每一个节点都是一个目标参数集合，若反演三个参数，则每个节点就是一个三行一列的矩阵；
%也可以将每个节点都理解为一个萤火虫，每只萤火虫都代表着一个目标反演参数数；再将萤火虫拉到同一行中，即每一列为一只萤火虫。
%接下来通过萤火虫算法让这些萤火虫进行发光和亮度发生变化（变亮可以理解为其反演的参数更优）。
for j=1:m
        for k=1:a
             for n=1:a^(m-j)
                 for i=1:a^(j-1) 
                 p(j,a^(j-1)*(k-1)+a^j*(n-1)+i) = number(j,k);
                 end
             end
         end
end  
N = size(p, 1);%参数数量 
M = size(p, 2);%列数即萤火虫数量

%Y1 = fModelCorrect + sqrt(0.01)*randn(I,1);  % measurement vector with guess noise
%Y2 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y3 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y = (Y1 + Y2 + Y3)./3;%多次测量取平均值以减少误差影响。

ITER=zeros(m,60);
ITER2=zeros(m,60);
ITER3=zeros(m,60);
d = 0;

ARF = 1.0;
gama = 1.0;
light0 = 1.0;
%BETA0 = 1.0;
% 数选取为目标似然函数的指数部分并用e来放大
%FBA = exp(-1*S)
% 现在定义萤火虫吸引力
%现在定义萤火虫移动规则及步长,每一列向量视作一个整体，

GaussDistri=sqrt(0.01)*randn(N,1);
%pbest = p();
k=0;
while (k<2)
    k=k+1;
    
    for i=1:M
            
      T = fModelGuess(p(:,i));
      T = T(1:NN);
      S = sum((Y-T).^2,1);
      S = sum(S,2);
      light = light0*gama*(1/S);                  %将适应度函数转化为光亮强度,gama为光吸收系数,用于扩大或者缩小适应度函数对光亮强度的影响
        if i == 1 
        lightmax = light;
        Imax = i;   
        lightmax2 = 0;
        Imax2 = i;  
        lightmax3 = 0;
        Imax3 = i;                               %,每一个参数的最大亮光参考值的标号
        else  
              if light>=lightmax
                lightmax = light;
                Imax = i;
              elseif light>=lightmax2
                lightmax2=light;
                Imax2 = i;
              elseif light>=lightmax3
                lightmax3=light;
                Imax3 = i;
              end      
        end
        
       
    end
    for i=1:M 
        r1=norm(p(j,Imax),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r2=norm(p(j,Imax2),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r3=norm(p(j,Imax3),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        BETA = lightmax*exp(-1*gama*(r1^2));%距离越远看到的亮度越小
        BETA2 = lightmax2*exp(-1*gama*(r2^2));%距离越远看到的亮度越小
        BETA3 = lightmax3*exp(-1*gama*(r3^2));%距离越远看到的亮度越小
        if BETA>=BETA2
            if BETA>=BETA3%一号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            else %三号最大
                pupdt(:,i)=p(:,i) + BETA3*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            end
        elseif BETA2>=BETA3%二号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax2)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        else  %三号大，此时二号第二
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        end
    end
    T=fModelGuess(p(:,i));
    T=T(1:NN);
    S = sum((Y-T).^2,1);
    S = sum(S,2);
    if S<1.0
        break
    end
    for i=1:m
    ITER(i,k)=p(i,Imax); 
    ITER2(i,k)=p(i,Imax2); 
    ITER3(i,k)=p(i,Imax3); 
    end 
   
end 


% Levenberg-Marquardt Method 
%pfirst = [p(1,Imax);p(2,Imax);p(3,Imax)];
p1 = p(:,Imax); % Guess parameters
p2 = p(:,Imax2); % Guess parameters
p3 = p(:,Imax3); % Guess parameters
p=p1;
pp=1;
N = size(p, 1);

%Y = fModelCorrect; 
%Y = load('celiang1.txt');

delp = 10^(-5);
J = zeros(I100, N);            % sensitivity matrix 
e1 = 3;
e2 = 1.5;
e3 = 3;
e4 = 5;
mu = 0.0001;

%plot(1:200, Y);

count = 0;
ite=0;
while true
 count = count + 1;
 if count>3
     break;
 end
 T = fModelGuess(p(:,i));
 T = T(1:NN);
 S = sum((Y-T).^2,1);
 S = sum(S,2);% Squared Error 

 %if norm(pupdt-p) < e4
         %delp = delp*10^(-1);
 %end
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         if abs(p(j)) < 0.001 %如果p值过小如0或者0.001，会导致雅可比矩阵为不可逆矩阵，从而导致迭代式子计算失败
         pInc(j) = pInc(j) + delp;
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         else
         pInc(j) = pInc(j) + delp*pInc(j);
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp*p(j));%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         end
         
     end
 end

 omega = diag(diag(J'*J));   

 pupdt = p;

 while true
     ite=ite+1;
     Yf=fModelGuess(pupdt);
     Yf=Yf(1:NN);
     pupdt = pupdt + (pinv(J'*J + mu*omega))*(J')*(Y - Yf);   %此处的阻尼因子mu*omega可优化
     T = fModelGuess(pupdt); 
     T = T(1:NN);        % updated T
     Supdt = sum((Y-T).^2);
     Supdt = sum(Supdt,2); 
     for i=1:m
     ITER(i,k+ite)=pupdt(i); 
     end
     if Supdt >= S         
       mu = 10*mu;
     elseif ite>=12
        p=p2;
        break;
     else
       mu = 0.1*mu;  
       break;
     end
    
 end
 T = fModelGuess(pupdt); 
 T = T(1:NN);
 if Supdt<e1
     pEstimate = pupdt;
     break;
 elseif ( norm(J'*(Y - T))<e2 || norm(pupdt-p)<e3)
    if pp==1
     p=p2;
     pp=pp+1;
    else
     p=p3;   
    end
 else
     p = pupdt;
 end
 
end

ps1 = pEstimate;






















NN = 200;%该部分代码的计算时间序列为0-100s
Y = Y200;
number = [-100 3 300;-10 5 50;-10 -1 0.1;];%待反演参数的取值范围
m = size(number,1);%待反演参数数量
a = size(number,2);
p = zeros(m,a^m);
%先创造参数空间，在这个参数空间里，每一个节点都是一个目标参数集合，若反演三个参数，则每个节点就是一个三行一列的矩阵；
%也可以将每个节点都理解为一个萤火虫，每只萤火虫都代表着一个目标反演参数数；再将萤火虫拉到同一行中，即每一列为一只萤火虫。
%接下来通过萤火虫算法让这些萤火虫进行发光和亮度发生变化（变亮可以理解为其反演的参数更优）。
for j=1:m
        for k=1:a
             for n=1:a^(m-j)
                 for i=1:a^(j-1) 
                 p(j,a^(j-1)*(k-1)+a^j*(n-1)+i) = number(j,k);
                 end
             end
         end
end  
N = size(p, 1);%参数数量 
M = size(p, 2);%列数即萤火虫数量

%Y1 = fModelCorrect + sqrt(0.01)*randn(I,1);  % measurement vector with guess noise
%Y2 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y3 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y = (Y1 + Y2 + Y3)./3;%多次测量取平均值以减少误差影响。

ITER=zeros(m,60);
ITER2=zeros(m,60);
ITER3=zeros(m,60);
d = 0;

ARF = 1.0;
gama = 1.0;
light0 = 1.0;
%BETA0 = 1.0;
% 数选取为目标似然函数的指数部分并用e来放大
%FBA = exp(-1*S)
% 现在定义萤火虫吸引力
%现在定义萤火虫移动规则及步长,每一列向量视作一个整体，

GaussDistri=sqrt(0.01)*randn(N,1);
%pbest = p();
k=0;
while (k<2)
    k=k+1;
    
    for i=1:M
            
      T = fModelGuess(p(:,i));
      T = T(1:NN);
      S = sum((Y-T).^2,1);
      S = sum(S,2);
      light = light0*gama*(1/S);                  %将适应度函数转化为光亮强度,gama为光吸收系数,用于扩大或者缩小适应度函数对光亮强度的影响
        if i == 1 
        lightmax = light;
        Imax = i;   
        lightmax2 = 0;
        Imax2 = i;  
        lightmax3 = 0;
        Imax3 = i;                               %,每一个参数的最大亮光参考值的标号
        else  
              if light>=lightmax
                lightmax = light;
                Imax = i;
              elseif light>=lightmax2
                lightmax2=light;
                Imax2 = i;
              elseif light>=lightmax3
                lightmax3=light;
                Imax3 = i;
              end      
        end
        
       
    end
    for i=1:M 
        r1=norm(p(j,Imax),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r2=norm(p(j,Imax2),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r3=norm(p(j,Imax3),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        BETA = lightmax*exp(-1*gama*(r1^2));%距离越远看到的亮度越小
        BETA2 = lightmax2*exp(-1*gama*(r2^2));%距离越远看到的亮度越小
        BETA3 = lightmax3*exp(-1*gama*(r3^2));%距离越远看到的亮度越小
        if BETA>=BETA2
            if BETA>=BETA3%一号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            else %三号最大
                pupdt(:,i)=p(:,i) + BETA3*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            end
        elseif BETA2>=BETA3%二号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax2)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        else  %三号大，此时二号第二
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        end
    end
    T=fModelGuess(p(:,i));
    T=T(1:NN);
    S = sum((Y-T).^2,1);
    S = sum(S,2);
    if S<1.0
        break
    end
    for i=1:m
    ITER(i,k)=p(i,Imax); 
    ITER2(i,k)=p(i,Imax2); 
    ITER3(i,k)=p(i,Imax3); 
    end 
   
end 


% Levenberg-Marquardt Method 
%pfirst = [p(1,Imax);p(2,Imax);p(3,Imax)];
p1 = p(:,Imax); % Guess parameters
p2 = p(:,Imax2); % Guess parameters
p3 = p(:,Imax3); % Guess parameters
p=p1;
pp=1;
N = size(p, 1);

%Y = fModelCorrect; 
%Y = load('celiang1.txt');

delp = 10^(-5);
J = zeros(I200, N);            % sensitivity matrix 
e1 = 3;
e2 = 1.5;
e3 = 3;
e4 = 5;
mu = 0.0001;

%plot(1:200, Y);

count = 0;
ite=0;
while true
 count = count + 1;
 if count>3
     break;
 end
 T = fModelGuess(p(:,i));
 T = T(1:NN);
 S = sum((Y-T).^2,1);
 S = sum(S,2);% Squared Error 

 %if norm(pupdt-p) < e4
         %delp = delp*10^(-1);
 %end
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         if abs(p(j)) < 0.001 %如果p值过小如0或者0.001，会导致雅可比矩阵为不可逆矩阵，从而导致迭代式子计算失败
         pInc(j) = pInc(j) + delp;
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         else
         pInc(j) = pInc(j) + delp*pInc(j);
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp*p(j));%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         end
         
     end
 end

 omega = diag(diag(J'*J));   

 pupdt = p;

 while true
     ite=ite+1;
     Yf=fModelGuess(pupdt);
     Yf=Yf(1:NN);
     pupdt = pupdt + (pinv(J'*J + mu*omega))*(J')*(Y - Yf);   %此处的阻尼因子mu*omega可优化
     T = fModelGuess(pupdt); 
     T = T(1:NN);        % updated T
     Supdt = sum((Y-T).^2);
     Supdt = sum(Supdt,2); 
     for i=1:m
     ITER(i,k+ite)=pupdt(i); 
     end
     if Supdt >= S         
       mu = 10*mu;
     elseif ite>=12
        p=p2;
        break;
     else
       mu = 0.1*mu;  
       break;
     end
    
 end
 T = fModelGuess(pupdt); 
 T = T(1:NN);
 if Supdt<e1
     pEstimate = pupdt;
     break;
 elseif ( norm(J'*(Y - T))<e2 || norm(pupdt-p)<e3)
    if pp==1
     p=p2;
     pp=pp+1;
    else
     p=p3;   
    end
 else
     p = pupdt;
 end
 
end

ps2 = pEstimate;




























NN = 300;%该部分代码的计算时间序列为0-100s
Y = Y300;
number = [-100 3 300;-10 5 50;-10 -1 0.1;];%待反演参数的取值范围
m = size(number,1);%待反演参数数量
a = size(number,2);
p = zeros(m,a^m);
%先创造参数空间，在这个参数空间里，每一个节点都是一个目标参数集合，若反演三个参数，则每个节点就是一个三行一列的矩阵；
%也可以将每个节点都理解为一个萤火虫，每只萤火虫都代表着一个目标反演参数数；再将萤火虫拉到同一行中，即每一列为一只萤火虫。
%接下来通过萤火虫算法让这些萤火虫进行发光和亮度发生变化（变亮可以理解为其反演的参数更优）。
for j=1:m
        for k=1:a
             for n=1:a^(m-j)
                 for i=1:a^(j-1) 
                 p(j,a^(j-1)*(k-1)+a^j*(n-1)+i) = number(j,k);
                 end
             end
         end
end  
N = size(p, 1);%参数数量 
M = size(p, 2);%列数即萤火虫数量

%Y1 = fModelCorrect + sqrt(0.01)*randn(I,1);  % measurement vector with guess noise
%Y2 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y3 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y = (Y1 + Y2 + Y3)./3;%多次测量取平均值以减少误差影响。

ITER=zeros(m,60);
ITER2=zeros(m,60);
ITER3=zeros(m,60);
d = 0;

ARF = 1.0;
gama = 1.0;
light0 = 1.0;
%BETA0 = 1.0;
% 数选取为目标似然函数的指数部分并用e来放大
%FBA = exp(-1*S)
% 现在定义萤火虫吸引力
%现在定义萤火虫移动规则及步长,每一列向量视作一个整体，

GaussDistri=sqrt(0.01)*randn(N,1);
%pbest = p();
k=0;
while (k<2)
    k=k+1;
    
    for i=1:M
            
      T = fModelGuess(p(:,i));
      T = T(1:NN);
      S = sum((Y-T).^2,1);
      S = sum(S,2);
      light = light0*gama*(1/S);                  %将适应度函数转化为光亮强度,gama为光吸收系数,用于扩大或者缩小适应度函数对光亮强度的影响
        if i == 1 
        lightmax = light;
        Imax = i;   
        lightmax2 = 0;
        Imax2 = i;  
        lightmax3 = 0;
        Imax3 = i;                               %,每一个参数的最大亮光参考值的标号
        else  
              if light>=lightmax
                lightmax = light;
                Imax = i;
              elseif light>=lightmax2
                lightmax2=light;
                Imax2 = i;
              elseif light>=lightmax3
                lightmax3=light;
                Imax3 = i;
              end      
        end
        
       
    end
    for i=1:M 
        r1=norm(p(j,Imax),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r2=norm(p(j,Imax2),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r3=norm(p(j,Imax3),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        BETA = lightmax*exp(-1*gama*(r1^2));%距离越远看到的亮度越小
        BETA2 = lightmax2*exp(-1*gama*(r2^2));%距离越远看到的亮度越小
        BETA3 = lightmax3*exp(-1*gama*(r3^2));%距离越远看到的亮度越小
        if BETA>=BETA2
            if BETA>=BETA3%一号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            else %三号最大
                pupdt(:,i)=p(:,i) + BETA3*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            end
        elseif BETA2>=BETA3%二号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax2)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        else  %三号大，此时二号第二
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        end
    end
    T=fModelGuess(p(:,i));
    T=T(1:NN);
    S = sum((Y-T).^2,1);
    S = sum(S,2);
    if S<1.0
        break
    end
    for i=1:m
    ITER(i,k)=p(i,Imax); 
    ITER2(i,k)=p(i,Imax2); 
    ITER3(i,k)=p(i,Imax3); 
    end 
   
end 


% Levenberg-Marquardt Method 
%pfirst = [p(1,Imax);p(2,Imax);p(3,Imax)];
p1 = p(:,Imax); % Guess parameters
p2 = p(:,Imax2); % Guess parameters
p3 = p(:,Imax3); % Guess parameters
p=p1;
pp=1;
N = size(p, 1);

%Y = fModelCorrect; 
%Y = load('celiang1.txt');

delp = 10^(-5);
J = zeros(I300, N);            % sensitivity matrix 
e1 = 3;
e2 = 1.5;
e3 = 3;
e4 = 5;
mu = 0.0001;

%plot(1:200, Y);

count = 0;
ite=0;
while true
 count = count + 1;
 if count>3
     break;
 end
 T = fModelGuess(p(:,i));
 T = T(1:NN);
 S = sum((Y-T).^2,1);
 S = sum(S,2);% Squared Error 

 %if norm(pupdt-p) < e4
         %delp = delp*10^(-1);
 %end
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         if abs(p(j)) < 0.001 %如果p值过小如0或者0.001，会导致雅可比矩阵为不可逆矩阵，从而导致迭代式子计算失败
         pInc(j) = pInc(j) + delp;
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         else
         pInc(j) = pInc(j) + delp*pInc(j);
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp*p(j));%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         end
         
     end
 end

 omega = diag(diag(J'*J));   

 pupdt = p;

 while true
     ite=ite+1;
     Yf=fModelGuess(pupdt);
     Yf=Yf(1:NN);
     pupdt = pupdt + (pinv(J'*J + mu*omega))*(J')*(Y - Yf);   %此处的阻尼因子mu*omega可优化
     T = fModelGuess(pupdt); 
     T = T(1:NN);        % updated T
     Supdt = sum((Y-T).^2);
     Supdt = sum(Supdt,2); 
     for i=1:m
     ITER(i,k+ite)=pupdt(i); 
     end
     if Supdt >= S         
       mu = 10*mu;
     elseif ite>=12
        p=p2;
        break;
     else
       mu = 0.1*mu;  
       break;
     end
    
 end
 T = fModelGuess(pupdt); 
 T = T(1:NN);

 if Supdt<e1
     pEstimate = pupdt;
     break;
 elseif ( norm(J'*(Y - T))<e2 || norm(pupdt-p)<e3)
    if pp==1
     p=p2;
     pp=pp+1;
    else
     p=p3;   
    end
 else
     p = pupdt;
 end
 
end

ps3 = pEstimate;
































NN = 400;%该部分代码的计算时间序列为0-100s
Y = Y400;
number = [-100 3 300;-10 5 50;-10 -1 0.1;];%待反演参数的取值范围
m = size(number,1);%待反演参数数量
a = size(number,2);
p = zeros(m,a^m);
%先创造参数空间，在这个参数空间里，每一个节点都是一个目标参数集合，若反演三个参数，则每个节点就是一个三行一列的矩阵；
%也可以将每个节点都理解为一个萤火虫，每只萤火虫都代表着一个目标反演参数数；再将萤火虫拉到同一行中，即每一列为一只萤火虫。
%接下来通过萤火虫算法让这些萤火虫进行发光和亮度发生变化（变亮可以理解为其反演的参数更优）。
for j=1:m
        for k=1:a
             for n=1:a^(m-j)
                 for i=1:a^(j-1) 
                 p(j,a^(j-1)*(k-1)+a^j*(n-1)+i) = number(j,k);
                 end
             end
         end
end  
N = size(p, 1);%参数数量 
M = size(p, 2);%列数即萤火虫数量

%Y1 = fModelCorrect + sqrt(0.01)*randn(I,1);  % measurement vector with guess noise
%Y2 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y3 = fModelCorrect + sqrt(0.01)*randn(I,1);          % measurement vector with guess noise
%Y = (Y1 + Y2 + Y3)./3;%多次测量取平均值以减少误差影响。

ITER=zeros(m,60);
ITER2=zeros(m,60);
ITER3=zeros(m,60);
d = 0;

ARF = 1.0;
gama = 1.0;
light0 = 1.0;
%BETA0 = 1.0;
% 数选取为目标似然函数的指数部分并用e来放大
%FBA = exp(-1*S)
% 现在定义萤火虫吸引力
%现在定义萤火虫移动规则及步长,每一列向量视作一个整体，

GaussDistri=sqrt(0.01)*randn(N,1);
%pbest = p();
k=0;
while (k<2)
    k=k+1;
    
    for i=1:M
            
      T = fModelGuess(p(:,i));
      T = T(1:NN);
      S = sum((Y-T).^2,1);
      S = sum(S,2);
      light = light0*gama*(1/S);                  %将适应度函数转化为光亮强度,gama为光吸收系数,用于扩大或者缩小适应度函数对光亮强度的影响
        if i == 1 
        lightmax = light;
        Imax = i;   
        lightmax2 = 0;
        Imax2 = i;  
        lightmax3 = 0;
        Imax3 = i;                               %,每一个参数的最大亮光参考值的标号
        else  
              if light>=lightmax
                lightmax = light;
                Imax = i;
              elseif light>=lightmax2
                lightmax2=light;
                Imax2 = i;
              elseif light>=lightmax3
                lightmax3=light;
                Imax3 = i;
              end      
        end
        
       
    end
    for i=1:M 
        r1=norm(p(j,Imax),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r2=norm(p(j,Imax2),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        r3=norm(p(j,Imax3),p(j,i));%计算最大亮度的向量和所有向量的范数，即两个向量间的距离，越大代表距离越远
                                     %如果自身够亮会不会也不会移动呢？即考虑了局部优解问题
        BETA = lightmax*exp(-1*gama*(r1^2));%距离越远看到的亮度越小
        BETA2 = lightmax2*exp(-1*gama*(r2^2));%距离越远看到的亮度越小
        BETA3 = lightmax3*exp(-1*gama*(r3^2));%距离越远看到的亮度越小
        if BETA>=BETA2
            if BETA>=BETA3%一号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            else %三号最大
                pupdt(:,i)=p(:,i) + BETA3*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
            end
        elseif BETA2>=BETA3%二号最大
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax2)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        else  %三号大，此时二号第二
                pupdt(:,i)=p(:,i) + BETA*(p(:,Imax3)-p(:,i));
                p(:,i) = pupdt(:,i);
               
        end
    end
    T=fModelGuess(p(:,i));
    T=T(1:NN);
    S = sum((Y-T).^2,1);
    S = sum(S,2);
    if S<1.0
        break
    end
    for i=1:m
    ITER(i,k)=p(i,Imax); 
    ITER2(i,k)=p(i,Imax2); 
    ITER3(i,k)=p(i,Imax3); 
    end 
   
end 


% Levenberg-Marquardt Method 
%pfirst = [p(1,Imax);p(2,Imax);p(3,Imax)];
p1 = p(:,Imax); % Guess parameters
p2 = p(:,Imax2); % Guess parameters
p3 = p(:,Imax3); % Guess parameters
p=p1;
pp=1;
N = size(p, 1);

%Y = fModelCorrect; 
%Y = load('celiang1.txt');

delp = 10^(-5);
J = zeros(I400, N);            % sensitivity matrix 
e1 = 3;
e2 = 1.5;
e3 = 3;
e4 = 5;
mu = 0.0001;

%plot(1:200, Y);

count = 0;
ite=0;
while true
 count = count + 1;
 if count>3
     break;
 end
 T = fModelGuess(p(:,i));
 T = T(1:NN);
 S = sum((Y-T).^2,1);
 S = sum(S,2);% Squared Error 

 %if norm(pupdt-p) < e4
         %delp = delp*10^(-1);
 %end
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         if abs(p(j)) < 0.001 %如果p值过小如0或者0.001，会导致雅可比矩阵为不可逆矩阵，从而导致迭代式子计算失败
         pInc(j) = pInc(j) + delp;
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         else
         pInc(j) = pInc(j) + delp*pInc(j);
         Tdelpj = fModelGuess(pInc);
         Tdelpj = Tdelpj(1:NN);
         Tpj = fModelGuess(p);
         Tpj =Tpj(1:NN);
         J(i,j) = (Tdelpj(i,1) - Tpj(i,1))/(delp*p(j));%计算雅可比矩阵，次数计算方式也许可进行优化
         %J(i+I,j) = (Tdelpj(i,2) - Tpj(i,2))/(delp);%计算雅可比矩阵，次数计算方式也许可进行优化
         end
         
     end
 end

 omega = diag(diag(J'*J));   

 pupdt = p;

 while true
     ite=ite+1;
     Yf=fModelGuess(pupdt);
     Yf=Yf(1:NN);
     pupdt = pupdt + (pinv(J'*J + mu*omega))*(J')*(Y - Yf);   %此处的阻尼因子mu*omega可优化
     T = fModelGuess(pupdt); 
     T = T(1:NN);        % updated T
     Supdt = sum((Y-T).^2);
     Supdt = sum(Supdt,2); 
     for i=1:m
     ITER(i,k+ite)=pupdt(i); 
     end
     if Supdt >= S         
       mu = 10*mu;
     elseif ite>=12
        p=p2;
        break;
     else
       mu = 0.1*mu;  
       break;
     end
    
 end
 T = fModelGuess(pupdt); 
 T = T(1:NN);
 if Supdt<e1
     pEstimate = pupdt;
     break;
 elseif ( norm(J'*(Y - T))<e2 || norm(pupdt-p)<e3)
    if pp==1
     p=p2;
     pp=pp+1;
    else
     p=p3;   
    end
 else
     p = pupdt;
 end
 
end

ps4 = pEstimate;














for i=1:400
 if i<101
  hflux(i)=ps1(1)+ps1(2)*i+ps1(3)*i*i;
 elseif i<201
    hflux(i)=ps2(1)+ps2(2)*i+ps2(3)*i*i;
 elseif i<301
    hflux(i)=ps3(1)+ps3(2)*i+ps3(3)*i*i;
 else
    hflux(i)=ps4(1)+ps4(2)*i+ps4(3)*i*i;
 end
end
t = 1:1:400;
yhflux=300*EXP(-(T-500)*(T-500)/1000)
plot(t,yhflux,'-r',t,hflux,'-xb');

plot(t,Y,'-r',t,T,'-xb');
Y1 = load('celiang2.txt');
T1 = fModelCorrect(p);
t = 1:1:400;
plot(t,Y1,'-r',t,T1,'-xb');

combos = cell(1,m);
[combos{:}]=ndgrid(number(:,1),number(:,2),number(:,3));
p = reshape(cat(m+1,combos{:}),m,[]);

%qt=485+14.*t+7.*t.^2;
%qtyuce=pEstimate(1)+t.*pEstimate(2)+ t.^2.*pEstimate(3);


% t1 = 4:0.1:5.5;
% qt1=30.*sin(t1).^1+0.*cos(t1)+0.*t1.^2;
% qtyuce1=cos(t1).*pEstimate(1)+sin(t1).*pEstimate(2)+ t1.^2.*pEstimate(3);

% figure;
% plot(t,qt,'-r',t,qtyuce,'-xb');
% xlabel('Time(s)');
% ylabel('Heat Flux(w/m2)');
% title('边界热流函数反演（σ²=0.1）')
% grid on;
% rectangle('Position', [min(t1) -30 max(t1)-min(t1) 10], 'EdgeColor','r');

% sub = axes('Position',[0.2,0.7,0.20,0.15]);
% plot(t1,qt1,'-r',t1,qtyuce1,'-xb');
% xlim([min(t1),max(t1)]);
% set(sub,'xtick',[],'ytick',[]);
