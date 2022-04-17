clear;
clear all;
clc;

%%%% ����������
lamda=632.8e-9; %%% �������Ⲩ��632.8nm
c=3e8;          %%% ����
k0=2*pi/lamda;  %%%%����650nm��k0Ϊ����
I0=10;          %%% �޹ⷴ������10mW
Ld=0.2;         %%% ������ǻ��0.2m
%%%% ��������
fs=1e5;             %%% ����Ƶ��
ts=1/fs;            %%% ��������
time=0:ts:0.03;      %%% ����ʱ��
%%%% �о��˶�����
L0=0.2;         %%% ��ʼ��ǻ��0.2m
d0=0.76e-6;        %%% �񶯷��ֵΪ2um
pha_amp=2*k0*d0; %%% ������λ���ֵΪ2um
f0=10;          %%% Ƶ��Ϊ10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------��ǻ���Ҳ���
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------��ǻ���ǲ���

% I_origin=I0*(1+m*cos(pha0));                                  %----------˫���������ź�


%%%% �ⷴ������
r4=0.21;
alfa=6;         %%% �߿�չ������
C=0.3;            %%% �ⷴ��ˮƽ����
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% ��ǿ����ϵ��


%%%% ����Ի�ϸ����źŹ�ǿ
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%��������
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------��ⷴ����λ
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------�Ի�ϸ����źŹ��ʽ�  
        end
    end
    index(i)=k;                                                 %----------��¼������λ��ĸ���
    i=i+1;
end

% %%%% ��ͼ������λ
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% ��ͼ��ǻԭʼ��λ   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% ��ͼ�Ի�ϸ����źŹ��ʽ�
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% ��ʼ���Ի�ϸ����ź�
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% ����ҽ�
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
lamda=632.8e-9; %%% �������Ⲩ��632.8nm
c=3e8;          %%% ����
k0=2*pi/lamda;  %%%%����650nm��k0Ϊ����
I0=10;          %%% �޹ⷴ������10mW
Ld=0.2;         %%% ������ǻ��0.2m
%%%% ��������
fs=1e5;             %%% ����Ƶ��
ts=1/fs;            %%% ��������
time=0:ts:0.03;      %%% ����ʱ��
%%%% �о��˶�����
L0=0.2;         %%% ��ʼ��ǻ��0.2m
d0=0.76e-6;        %%% �񶯷��ֵΪ2um
pha_amp=2*k0*d0; %%% ������λ���ֵΪ2um
f0=10;          %%% Ƶ��Ϊ10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------��ǻ���Ҳ���
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------��ǻ���ǲ���

% I_origin=I0*(1+m*cos(pha0));                                  %----------˫���������ź�


%%%% �ⷴ������
r4=0.21;
alfa=6;         %%% �߿�չ������
C=0.9;            %%% �ⷴ��ˮƽ����
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% ��ǿ����ϵ��


%%%% ����Ի�ϸ����źŹ�ǿ
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%��������
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------��ⷴ����λ
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------�Ի�ϸ����źŹ��ʽ�  
        end
    end
    index(i)=k;                                                 %----------��¼������λ��ĸ���
    i=i+1;
end

% %%%% ��ͼ������λ
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% ��ͼ��ǻԭʼ��λ   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% ��ͼ�Ի�ϸ����źŹ��ʽ�
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% ��ʼ���Ի�ϸ����ź�
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% ����ҽ�
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
lamda=632.8e-9; %%% �������Ⲩ��632.8nm
c=3e8;          %%% ����
k0=2*pi/lamda;  %%%%����650nm��k0Ϊ����
I0=10;          %%% �޹ⷴ������10mW
Ld=0.2;         %%% ������ǻ��0.2m
%%%% ��������
fs=1e5;             %%% ����Ƶ��
ts=1/fs;            %%% ��������
time=0:ts:0.03;      %%% ����ʱ��
%%%% �о��˶�����
L0=0.2;         %%% ��ʼ��ǻ��0.2m
d0=0.76e-6;        %%% �񶯷��ֵΪ2um
pha_amp=2*k0*d0; %%% ������λ���ֵΪ2um
f0=10;          %%% Ƶ��Ϊ10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------��ǻ���Ҳ���
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------��ǻ���ǲ���

% I_origin=I0*(1+m*cos(pha0));                                  %----------˫���������ź�


%%%% �ⷴ������
r4=0.21;
alfa=6;         %%% �߿�չ������
C=1.5;            %%% �ⷴ��ˮƽ����
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% ��ǿ����ϵ��


%%%% ����Ի�ϸ����źŹ�ǿ
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%��������
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------��ⷴ����λ
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------�Ի�ϸ����źŹ��ʽ�  
        end
    end
    index(i)=k;                                                 %----------��¼������λ��ĸ���
    i=i+1;
end

% %%%% ��ͼ������λ
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% ��ͼ��ǻԭʼ��λ   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% ��ͼ�Ի�ϸ����źŹ��ʽ�
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% ��ʼ���Ի�ϸ����ź�
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% ����ҽ�
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
lamda=632.8e-9; %%% �������Ⲩ��632.8nm
c=3e8;          %%% ����
k0=2*pi/lamda;  %%%%����650nm��k0Ϊ����
I0=10;          %%% �޹ⷴ������10mW
Ld=0.2;         %%% ������ǻ��0.2m
%%%% ��������
fs=1e5;             %%% ����Ƶ��
ts=1/fs;            %%% ��������
time=0:ts:0.03;      %%% ����ʱ��
%%%% �о��˶�����
L0=0.2;         %%% ��ʼ��ǻ��0.2m
d0=0.76e-6;        %%% �񶯷��ֵΪ2um
pha_amp=2*k0*d0; %%% ������λ���ֵΪ2um
f0=10;          %%% Ƶ��Ϊ10Hz
L=d0*2*pi*f0*(time+0.2);% L=d0*sin(2*pi*f0*time-pi/2+0.2);
pha0=2*k0*(L0+L);                            %----------��ǻ���Ҳ���
% pha0=2*k0*(L0+d0*sawtooth(2*pi*f0*time,0.5));                 %----------��ǻ���ǲ���

% I_origin=I0*(1+m*cos(pha0));                                  %----------˫���������ź�


%%%% �ⷴ������
r4=0.21;
alfa=6;         %%% �߿�չ������
C=2.1;            %%% �ⷴ��ˮƽ����
% C=L0/Ld*sqrt(1+alfa^2)*r4*0.81;
m=0.01;         %%% ��ǿ����ϵ��


%%%% ����Ի�ϸ����źŹ�ǿ
i=1;
for t=time
    n_temp=[];ph_temp=[];error_temp=[];%%%%��������
    [n_temp,pha_temp,error_temp]=solve_pha(pha0(i),C,alfa);     %----------��ⷴ����λ
    
    k=0;
    for j=1:n_temp
        if abs(error_temp(j))<=1e-3
            k=k+1;
            pha(i,k)=pha_temp(j);
            I_feedback(i,k)=I0+m*I0*cos(pha(i,k));              %----------�Ի�ϸ����źŹ��ʽ�  
        end
    end
    index(i)=k;                                                 %----------��¼������λ��ĸ���
    i=i+1;
end

% %%%% ��ͼ������λ
% subplot(211);       
% plot(time,pha0);
% for i=1:length(time)
%    for j=1:index(i)
%        hold on;plot(time(i),pha(i,j));
%    end
% end
%    
% %%%% ��ͼ��ǻԭʼ��λ   
% subplot(212);
% plot(time,I_origin);
% figure;

%%%% ��ͼ�Ի�ϸ����źŹ��ʽ�
% for i=1:length(time)
%     j=1:index(i);
%     hold on;plot(i*ones(1,index(i)),I_feedback(i,1:index(i)),'r*');
% end

OOP=zeros(1,length(time));            %%% ��ʼ���Ի�ϸ����ź�
OOP(1)=I_feedback(1,1);
for i=2:length(time)                  %%% ����ҽ�
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