addpath('quaternion_library');      % include quaternion library
close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal
%% Import and plot sensor data
Data = importdata('NAV_2.mat');
gx=Data(:,27); %%陀螺仪数据
gy=Data(:,28); 
gz=Data(:,29); 
ax=Data(:,9); %%加速度计数据
ay=Data(:,10); 
az=Data(:,11); 
mx=Data(:,15); %%磁强计数据
my=Data(:,16); 
mz=Data(:,17); 
rad2deg=180/pi;
deg2rad=pi/180;
time(1)=0;
quaternion = zeros(length(time), 4);
%% Process sensor data through algorithm
Pitch0=asin(ax(1)/sqrt(ax(1)*ax(1)+ay(1)*ay(1)+az(1)*az(1)));%初始角计算
Roll0=atan2(-ay(1),-az(1));
Yaw0=atan2(-my(1)*cos(Roll0)+mz(1)*sin(Roll0),mx(1)*cos(Pitch0)+my(1)*sin(Pitch0)*sin(Roll0)+mz(1)*sin(Pitch0)*cos(Roll0))-8.3*pi/180;
[q0,q1,q2,q3]=Quaternion_FromEulerAngle(Roll0,Pitch0,Yaw0);%初始四元数计算
q(1,:)=[q0,q1,q2,q3];
Pitch(1)= Pitch0*rad2deg; 
Roll(1) =Roll0*rad2deg;
Yaw(1) = Yaw0*rad2deg;
tt(1)=0;                                      %算法单次仿真时间测量
SamplePeriod=0.02;                            %数据采样时间
Beta=0.009;                                    %梯度下降法梯度步长设置
Kp=2;Ki=0.1;                                  %叉积法比例积分系数设置
eInt(1,:)=[0 0 0];                            %叉积法误差初值
Length= size(Data,1);                        
Re=6378137;r=6356752.3142;f=1/298.257;e=0.0818;wie=7.292e-5;%Re长半轴 r短半轴 f椭球扁率 e椭球偏心率 wie地球自转角速率
Ven(:,1)=[Data(1,6) Data(1,7) Data(1,8)]';
E(1)=Data(1,2)*deg2rad;L(1)=Data(1,3)*deg2rad;H(1)=Data(1,4);
fn(:,1)=[0 0 0]';g(:,1)=[0 0 0]';Venq1(:,1)=[0 0 0]'; Venq2(:,1)=[0 0 0]';
g0=9.7803;
%惯导相关的噪声统计数据
Q_wg  = (1/(57*3600))^2;         %陀螺的随机漂移为0.5度每小时
Q_wa  = ((0.5e-4)*g0)^2;                  %加速度计的随机偏差为0.5e-4*g  
Q 		= diag([Q_wg Q_wg Q_wg,  Q_wa Q_wa Q_wa]);%系统噪声方差
%Q=diag([1e-9 1e-9 1e-9 1e-7 1e-7 1e-7 ]);

Tg 		= 300*ones(3,1);                          %Tg     陀螺仪误差漂移相关时间
Ta 		= 1000*ones(3,1);                         %Ta     加表误差漂移相关时间
%error_cz=0.01*pi/180/3600;%%%陀螺相关漂移 0.5度/时，
%tg=300; %%%相关时间tg=300s
%kt=sqrt(2*error_cz^2/tg);%%%驱动噪声强度均方根
%ta=1000;%%相关时间
%Ra=1.0e-4*g0;%%均方根Ra=0.0001*g
%ka=sqrt(2*Ra^2/ta);%%%%驱动噪声强度均方根
%zt=0.05*deg2rad*ones(1,3);  %%%%零均值姿态误差噪声均方根向量 
%wz=[0.001/60*pi/180 0.001/60*pi/180 0.1];   %%%%零均值位置噪声均方根向量 
%sd=0.0001*[1,1,1];    %%%%%零均值速度噪声均方根向量 
%fc=[zt sd wz kt kt kt ka ka ka ].^2;%%%%%新的方差阵对角线元素
%Q=diag(fc);
%R = 1e-6*eye(6);                               %GPS误差，测量噪声方差
Rlamt=1e-5*pi/(60*180); %%%经纬度误差均方根,弧度制
Rl=1e-5*pi/(60*180);
Rh=1e-11; %%%高度误差均方根,单位 米
Rvx=1e-7; %%%速度误差均方根,单位 米/秒
Rvy=1e-7;
Rvz=5e-9;
K=[Rlamt Rl Rh Rvx Rvy Rvz];%
R=diag(K);
m=0;
PP0(1:15,1:15,1) = diag([0.1/(57) 0.1/(57) 0.1/57, 0.01 0.01 0.01, 0 0 0, 0.1/(57*3600) 0.1/(57*3600) 0.1/(57*3600), (1e-4)*g0 (1e-4)*g0 (1e-4)*g0].^2);   %初始误差协方差阵,加速度计的初始偏值均取1e-4*g 陀螺的常值漂移取0.1度每小时
X(:,:,1)= zeros(15,1);                            %初始状态
 a=0;                               %tao    迭代步长
E_pv2(1,:)=[0 0 0 0 0 0];
for t = 2:Length
    tic;
    q(t,:) = madgwickAHRS(Data(t,27:29), Data(t,9:11), Data(t,15:17),q(t-1, :),Beta,SamplePeriod);	%梯度下降法  
   % [q(t,:),eInt(t,:)] = MahonyAHRSupdate( Data(t,27:29), Data(t,9:11), Data(t,15:17),q(t-1, :),SamplePeriod,Ki,Kp,eInt(t-1,:));%Mahony法
    tt(t)=toc;%计算算法运行时间
    
    T=[ 1 - 2 * (q(t,4) *q(t,4) + q(t,3) * q(t,3)) 2 * (q(t,2) * q(t,3) +q(t,1) * q(t,4)) 2 * (q(t,2) * q(t,4)-q(t,1) * q(t,3));
        2 * (q(t,2) * q(t,3)-q(t,1) * q(t,4)) 1 - 2 * (q(t,4) *q(t,4) + q(t,2) * q(t,2)) 2 * (q(t,3) * q(t,4)+q(t,1) * q(t,2));
        2 * (q(t,2) * q(t,4) +q(t,1) * q(t,3)) 2 * (q(t,3) * q(t,4)-q(t,1) * q(t,2)) 1 - 2 * (q(t,2) *q(t,2) + q(t,3) * q(t,3))];%cnb
    Roll(t)  = atan2(T(2,3),T(3,3))*rad2deg;
    Pitch(t) = asin(-T(1,3))*rad2deg;
    Yaw(t)   = atan2(T(1,2),T(1,1))*rad2deg-8.3;  
    
    %Rm(t-1)=Re*(1-e*e)/power((1-e*e*sin(L(t-1))*sin(L(t-1))),1.5);%子午曲率半径
   % Rn(t-1)=Re/sqrt(1-e*e*sin(L(t-1))*sin(L(t-1)));               %卯酉曲率半径
    Rm(t-1)=Re*(1-2*f+3*f*sin(L(t-1))*sin(L(t-1))); 
    Rn(t-1)=Re*(1+f*sin(L(t-1))*sin(L(t-1)));
    R0=sqrt(Rm(t-1)*Rn(t-1));                                     %曲率平均半径
    Wien(t-1,:)=[wie*cos(L(t-1)) 0 -wie*sin(L(t-1))];             %北东地坐标系中地球自转角速度
    Wenn(t-1,:)=[Ven(2,t-1)/(Rn(t-1)+H(t-1)) -Ven(1,t-1)/(Rm(t-1)+H(t-1)) -Ven(2,t-1)*tan(L(t-1))/(Rm(t-1)+H(t-1))];%地理坐标系相对于地球固连坐标系的转动角速率
    g0=9.780318;%*(1+5.3024*1e-3*sin(L(t-1))*sin(L(t-1))-5.9*1e-6*sin(2*L(t-1))*sin(2*L(t-1)));
    %g(:,t)=[0 0 g0/(1+H(t-1)/R0)^2]';
    %g(:,t)=[0 0 g0*(1-2*H(t-1)/Re)]';
    Fn(:,t-1)=T'*([Data(t-1,9:11)]')*9.8;                         %分解到地理坐标系中的比力
    WWX=[0 2*Wien(t-1,3)+Wenn(t-1,3) -(2*Wien(t-1,2)+Wenn(t-1,2));-(2*Wien(t-1,3)+Wenn(t-1,3)) 0 2*Wien(t-1,1)+Wenn(t-1,1);2*Wien(t-1,2)+Wenn(t-1,2) -(2*Wien(t-1,1)+Wenn(t-1,1)) 0];%比力误差项
    Venq1(:,t)=Fn(:,t-1)+WWX*Ven(:,t-1)+[0 0 10.05]';                    %速度微分
    Venq2(:,t)=Data(t-1,24:26)';
   
    Ven(:,t)=Ven(:,t-1)+Venq1(:,t)*SamplePeriod;    
    L(t)=L(t-1)+(Ven(1,t-1)/(Rm(t-1)+H(t-1)))*SamplePeriod;              %纬度
    E(t)=E(t-1)+(Ven(2,t-1)/(cos(L(t-1))*(Rn(t-1)+H(t-1))))*SamplePeriod;%经度
    H(t)=H(t-1)-Ven(3,t-1)*SamplePeriod;                                 %高度
    time(t)=time(t-1)+0.02;
    
    m=m+1;
    if Data(t,36)==1
        a=a+1;
        tao=m*0.02;
        m=0;
        Dpv=[L(t)*rad2deg-Data(t,3),E(t)*rad2deg-Data(t,2),H(t)-Data(t,4),Ven(1,t)-Data(t,6),Ven(2,t)-Data(t,7),Ven(3,t)-Data(t,8)];       %惯导和Gps所测位置速度之差
        [E_v, E_p,PP,XX] = kalman_GPS_INS_pv(Dpv, Ven(:,t), L(t),H(t), T, Fn(:,t-1), Q, R, Tg, Ta, tao,Rm(t-1),Rn(t-1),PP0(:,:,a),X(:,:,a)) ;%GPS/INS位置速度组合 卡尔曼滤波
        PP0(:,:,a+1)=PP;
        X(:,:,a+1)=XX;
        %{
        if t~=Length
        Data(t+1,9:11)=Data(t+1,9:11)-X(13:15,a+1)';
        Data(t+1,27:29)=Data(t+1,27:29)-X(10:12,a+1)';
        end
        %}
        % Q=Q1;
        Ven(:,t)=Ven(:,t)-[E_v(1); E_v(2); 0.2*E_v(3)];
        L(t)=L(t)-0.29*E_p(1);
        E(t)=E(t)-0.32*E_p(2);
        H(t)=H(t)-E_p(3);
   end   
  % E_pv2(t,:)=[Ven(1,t)-Data(t,21),Ven(2,t)-Data(t,22),Ven(3,t)-Data(t,23),E(t)*rad2deg-Data(t,18),L(t)*rad2deg-Data(t,19),H(t)-Data(t,20)]; 
end
%save E_pv2;
%save Ven;save L;save E;save H;
for i=1:Length
   L(i)=L(i)*rad2deg; 
    E(i)=E(i)*rad2deg; 
end
figure('Name', 'Euler Angles');
hold on;
plot(time, Roll, 'r',time,Data(:,30)*57.3,'m');
plot(time, Pitch, 'g',time,Data(:,31)*57.3,'c');
plot(time, Yaw, 'b',time,Data(:,32)*57.3,'k');
title('Euler angles');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('Roll计算值','Roll参考值', 'Pitch计算值', 'Pitch参考值','Yaw计算值','Yaw参考值');
hold off;
figure('Name', '1');
hold on;
plot(time,E, 'r');
plot(time, Data(:,18), 'g');
title('经度对比图','fontsize',20);
xlabel('Time (s)','fontsize',16);
ylabel(' longitude','fontsize',16);
d=legend('经度计算值', '经度参考值');
set(d,'Fontsize',10,'LineWidth',2,'box','off');
hold off;
figure('Name', '2');
hold on;
plot(time, L, 'r');
plot(time, Data(:,19), 'g');
title('纬度对比图','fontsize',20);
xlabel('Time (s)','fontsize',16);
ylabel('latitude','fontsize',16);
d=legend('纬度计算值', '纬度参考值');
set(d,'Fontsize',10,'LineWidth',2,'box','off');
hold off;
figure('Name', '3');
hold on;
plot(time, H, 'r');
plot(time, Data(:,20), 'g');
title('高度对比图','fontsize',20);
xlabel('Time (s)','fontsize',16);
ylabel(' H','fontsize',16);
d=legend('高度计算值', '高度参考值');
set(d,'Fontsize',10,'LineWidth',2,'box','off');
hold off;
figure('Name', '4');
hold on;
plot(time, Ven(1,:), 'r');
plot(time, Data(:,21), 'g');
title('北向速度对比图','fontsize',20);
xlabel('Time (s)','fontsize',16);
ylabel('北向速度','fontsize',16);
d=legend('北向速度计算值', '北向速度参考值');
set(d,'Fontsize',10,'LineWidth',2,'box','off');
hold off;
figure('Name', '5');
hold on;
plot(time,Ven(2,:), 'r');
plot(time, Data(:,22), 'g');
title('东向速度对比图','fontsize',20);
xlabel('Time (s)','fontsize',16);
ylabel('东向速度','fontsize',16);
d=legend('东向速度计算值', '东向速度参考值');
set(d,'Fontsize',10,'LineWidth',2,'box','off');
hold off;
figure('Name', '6');
hold on;
plot(time, Ven(3,:), 'r');
plot(time, Data(:,23), 'g');
title('地向速度对比图','fontsize',20);
xlabel('Time(s)','fontsize',16);
ylabel('地向速度','fontsize',16);
d=legend('地向速度计算值', '地向速度参考值');
set(d,'Fontsize',10,'LineWidth',2,'box','off');
hold off;
