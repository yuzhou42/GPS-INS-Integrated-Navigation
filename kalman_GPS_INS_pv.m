%GPS/INS位置速度组合 卡尔曼滤波
%连续状态系统方程
%dx = F*x + G*w             %z = H*x + v
%离散状态系统方程
%x(k+1) = A*x(k) + B*w(k)   %z(k+1) = C*x(k+1) + v(k+1)
%Dpv     量测位置速度误差， 作为滤波器输入
%Q      系统噪声方差
%R      测量噪声方差
%Ta     加表误差漂移相关时间
%Tg     陀螺仪误差漂移相关时间
%tao    迭代步长
function [E_v, E_p, PP,XX] = kalman_GPS_INS_pv(Dpv, Ven, L,h, T, Fn, Q, R, Tg, Ta, tao,Rm,Rn,PP0,X)
    wie 	= 7.292e-5;  %地球自转角速度
    V=Ven;
    Vn 		= V(1);
    Ve 		= V(2);
    Vd		= V(3);
    Fn=Fn;
    fn 		= Fn(1);
    fe 		= Fn(2);
    fd 		= Fn(3);
    %连续系统状态转换阵 F 的时间更新
    F = zeros(15,15);
    F(1,2)       = -wie*sin(L)-Ve*tan(L)/(Rn+h);
    F(1,3)       = Vn/(Rm+h);
    F(1,5)       = 1/(Rn+h);
    F(1,7)       = -wie*sin(L);
    F(1,9)       = -Ve/(Rn+h)^2;
    F(2,1)       = wie*sin(L)+Ve*tan(L)/(Rn+h);
    F(2,3)       = wie*cos(L)+Ve/(Rn+h);
    F(2,4)       = -1/(Rm+h);
    F(2,9)       = Vn/(Rm+h)^2;
    F(3,1)       = -Vn/(Rm+h);
    F(3,2)       = -wie*cos(L)-Ve/(Rn+h);
    F(3,5)       =-tan(L)/(Rn+h);
    F(3,7)       = -wie*cos(L)-Ve*(sec(L)^2)/(Rn+h);
    F(3,9)       = Ve*tan(L)/(Rn+h)^2;
    F(4,2)       = -fd;
    F(4,3)       = fe;
    F(4,4)       =Vd/(Rm+h);
    F(4,5)       =-2*wie*sin(L)-Ve*tan(L)/(Rn+h)-Ve*tan(L)/(Rn+h);
    F(4,6)       = Vn/(Rm+h);
    F(4,7)       = -2*wie*cos(L)*Ve-Ve^2*sec(L)^2/(Rn+h);
    F(4,9)       = Ve^2*tan(L)/(Rn+h)^2-Vn*Vd/(Rm+h)^2;
    F(5,1)       = fd;
    F(5,3)       = -fn;
    F(5,4)       = 2*wie*sin(L)+Ve*tan(L)/(Rn+h);
    F(5,5)       = (Vn*tan(L)+Vd)/(Rn+h);
    F(5,6)       = 2*wie*cos(L)+Ve/(Rn+h);
    F(5,7)       = (2*wie*cos(L)+Ve*(sec(L)^2)/(Rn+h))*Vn-2*wie*Vd*sin(L);
    F(5,9)       = (Ve*Vn*tan(L)-Ve*Vd)/(Rn+h)^2;
    F(6,1)       = -fe;
    F(6,2)       = fn;
    F(6,4)       = -Vn/(Rm+h);
    F(6,5)       = -2*wie*cos(L)-2*Ve/(Rn+h);
    F(6,7)       = 2*Ve*wie*sin(L);
    F(6,9)       = Ve^2/(Rn+h)^2+Vn^2/(Rm+h)^2;
    F(7,4)       = 1/(Rm+h);
    F(7,9)       = -Vn/(Rm+h)^2;
    F(8,5)       = 1/((Rn+h)*cos(L));
    F(8,7)       = Ve*tan(L)/((Rn+h)*cos(L));
    F(8,9)       = -Ve/(cos(L)*(Rn+h)^2);
    F(9,6)       = -1;
    F(1:3,10:12) = T';
    F(4:6,13:15) = T';
    
    F(10,10)     = -1/Tg(1);
    F(11,11)     = -1/Tg(2);
    F(12,12)     = -1/Tg(3);
    F(13,13)     = -1/Ta(1);
    F(14,14)     = -1/Ta(2);
    F(15,15)     = -1/Ta(3);
    
    %连续系统输入矩阵更新
    G = zeros(15,6);
    G(1:3,1:3)   = T';
    G(4:6,4:6) =  T';  
    %连续系统量测阵更新
    H = zeros(6,15);
    H(1,7) 			= Rm+h;
    H(2,8) 			= (Rn+h)*cos(L);
    H(3,9) 			= 1;
    H(4,4) 			= 1;
    H(5,5) 			= 1;
    H(6,6) 			= 1;
    %连续系统离散化
    A= eye(15,15)+F*tao;
    %B= (eye(15,15)+tao*F/2)*G*tao; 
    B= G*tao;
    %卡尔曼滤波
    X=A*X;
    P= A*(PP0)*A'+B*Q*B';
    %P= A*(PP0)*A'+Q;
    K= P*H'*inv(H*P*H'+R);
    %PP= (eye(15,15)-K*H)*P*(eye(15,15)-K*H)'+K*R*K';
    PP= (eye(15,15)-K*H)*P;
    z = Dpv';
    XX = X+K*(z-H*X);
    %Q1=(Q+A*Q*A')*tao/2;    %%循环处理
    E_v= XX(4:6);
    E_p= XX(7:9);

