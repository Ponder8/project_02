clear;
clear all;
clc;

%%%% 激光器参数
lamda=632.8e-9; %%% 激光器光波长632.8nm
c=3e8;          %%% 光速
k0=2*pi/lamda;  %%%%波长650nm，k0为波数
I0=10;          %%% 无光反馈功率10mW
Ld=0.2;         %%% 激光器腔长0.2m
%%%% 采样参数
fs=1e5;             %%% 采样频率
ts=1/fs;            %%% 采样周期
time=0:ts:0.03;      %%% 测量时间
%%%% 靶镜运动参数
L0=0.2;         %%% 初始外腔长0.2m
d0=0.76e-6;        %%% 振动峰峰值为2um
pha_amp=2*k0*d0; %%% 正弦相位峰峰值为2um
f0=10;          %%% 频率为10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------外腔正弦波振动
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------外腔三角波振动

% I_origin=I0*(1+m*cos(pha0));                                  %----------双光束干涉信号


%%%% 光反馈参数
r4=0.21;
alfa=6;         %%% 线宽展宽因子
C=0.3;            %%% 光反馈水平因子
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% 光强调制系数


%%%% 求解自混合干涉信号光强
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%清空求解结果
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------求解反馈相位
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------自混合干涉信号功率解  
        end
    end
    index(i)=k;                                                 %----------记录反馈相位解的个数
    i=i+1;
end

% %%%% 画图反馈相位
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% 画图外腔原始相位   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% 画图自混合干涉信号功率解
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% 初始化自混合干涉信号
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% 向后找解
    dev_temp=[];
    for j=1:index(i)
       dev_temp(j)=abs(I_feedback(i,j)-OOP(i-1)); 
    end
    [min_value,min_index]=min(dev_temp);
    OOP(i)=I_feedback(i,min_index);
end


zz=time/20*1e4;
qq=(OOP-mean(OOP))*50/6;
figure,plot(zz,qq,'r');
hold on;
lamda=632.8e-9; %%% 激光器光波长632.8nm
c=3e8;          %%% 光速
k0=2*pi/lamda;  %%%%波长650nm，k0为波数
I0=10;          %%% 无光反馈功率10mW
Ld=0.2;         %%% 激光器腔长0.2m
%%%% 采样参数
fs=1e5;             %%% 采样频率
ts=1/fs;            %%% 采样周期
time=0:ts:0.03;      %%% 测量时间
%%%% 靶镜运动参数
L0=0.2;         %%% 初始外腔长0.2m
d0=0.76e-6;        %%% 振动峰峰值为2um
pha_amp=2*k0*d0; %%% 正弦相位峰峰值为2um
f0=10;          %%% 频率为10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------外腔正弦波振动
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------外腔三角波振动

% I_origin=I0*(1+m*cos(pha0));                                  %----------双光束干涉信号


%%%% 光反馈参数
r4=0.21;
alfa=6;         %%% 线宽展宽因子
C=0.9;            %%% 光反馈水平因子
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% 光强调制系数


%%%% 求解自混合干涉信号光强
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%清空求解结果
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------求解反馈相位
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------自混合干涉信号功率解  
        end
    end
    index(i)=k;                                                 %----------记录反馈相位解的个数
    i=i+1;
end

% %%%% 画图反馈相位
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% 画图外腔原始相位   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% 画图自混合干涉信号功率解
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% 初始化自混合干涉信号
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% 向后找解
    dev_temp=[];
    for j=1:index(i)
       dev_temp(j)=abs(I_feedback(i,j)-OOP(i-1)); 
    end
    [min_value,min_index]=min(dev_temp);
    OOP(i)=I_feedback(i,min_index);
end


zz=time/20*1e4;
ee=(OOP-mean(OOP))*50/6;
plot(zz,ee,'b');
hold on
lamda=632.8e-9; %%% 激光器光波长632.8nm
c=3e8;          %%% 光速
k0=2*pi/lamda;  %%%%波长650nm，k0为波数
I0=10;          %%% 无光反馈功率10mW
Ld=0.2;         %%% 激光器腔长0.2m
%%%% 采样参数
fs=1e5;             %%% 采样频率
ts=1/fs;            %%% 采样周期
time=0:ts:0.03;      %%% 测量时间
%%%% 靶镜运动参数
L0=0.2;         %%% 初始外腔长0.2m
d0=0.76e-6;        %%% 振动峰峰值为2um
pha_amp=2*k0*d0; %%% 正弦相位峰峰值为2um
f0=10;          %%% 频率为10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------外腔正弦波振动
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------外腔三角波振动

% I_origin=I0*(1+m*cos(pha0));                                  %----------双光束干涉信号


%%%% 光反馈参数
r4=0.21;
alfa=6;         %%% 线宽展宽因子
C=1.5;            %%% 光反馈水平因子
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% 光强调制系数


%%%% 求解自混合干涉信号光强
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%清空求解结果
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------求解反馈相位
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------自混合干涉信号功率解  
        end
    end
    index(i)=k;                                                 %----------记录反馈相位解的个数
    i=i+1;
end

% %%%% 画图反馈相位
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% 画图外腔原始相位   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% 画图自混合干涉信号功率解
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% 初始化自混合干涉信号
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% 向后找解
    dev_temp=[];
    for j=1:index(i)
       dev_temp(j)=abs(I_feedback(i,j)-OOP(i-1)); 
    end
    [min_value,min_index]=min(dev_temp);
    OOP(i)=I_feedback(i,min_index);
end


zz=time/20*1e4;

ss=(OOP-mean(OOP))*50/6;
plot(zz,ss,'g');
hold on;
lamda=632.8e-9; %%% 激光器光波长632.8nm
c=3e8;          %%% 光速
k0=2*pi/lamda;  %%%%波长650nm，k0为波数
I0=10;          %%% 无光反馈功率10mW
Ld=0.2;         %%% 激光器腔长0.2m
%%%% 采样参数
fs=1e5;             %%% 采样频率
ts=1/fs;            %%% 采样周期
time=0:ts:0.03;      %%% 测量时间
%%%% 靶镜运动参数
L0=0.2;         %%% 初始外腔长0.2m
d0=0.76e-6;        %%% 振动峰峰值为2um
pha_amp=2*k0*d0; %%% 正弦相位峰峰值为2um
f0=10;          %%% 频率为10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------外腔正弦波振动
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------外腔三角波振动

% I_origin=I0*(1+m*cos(pha0));                                  %----------双光束干涉信号


%%%% 光反馈参数
r4=0.21;
alfa=6;         %%% 线宽展宽因子
C=2.1;            %%% 光反馈水平因子
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% 光强调制系数


%%%% 求解自混合干涉信号光强
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%清空求解结果
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------求解反馈相位
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------自混合干涉信号功率解  
        end
    end
    index(i)=k;                                                 %----------记录反馈相位解的个数
    i=i+1;
end

% %%%% 画图反馈相位
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% 画图外腔原始相位   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% 画图自混合干涉信号功率解
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% 初始化自混合干涉信号
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% 向后找解
    dev_temp=[];
    for j=1:index(i)
       dev_temp(j)=abs(I_feedback(i,j)-OOP(i-1)); 
    end
    [min_value,min_index]=min(dev_temp);
    OOP(i)=I_feedback(i,min_index);
end


zz=time/20*1e4;
qq=(OOP-mean(OOP))*50/6;
plot(zz,qq,'k');