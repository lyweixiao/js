function [x,y]=DrawPoint(x0,y0,R,num_Dian);
% ���������ǣ�Բ�ĺ������꣬�뾶�͵������
theta=0:0.001:360;
% ���ü�����õ�Բ������
Circle1=x0+R*cos(theta);
Circle2=y0+R*sin(theta);
% ��Բ
% plot(Circle1,Circle2,'r');
hold on;
% �������num_Dian���뾶
r=R*sqrt(rand(1,num_Dian));
% �õ����ɵ�ĽǶȣ������ü�������ʽ������
seta=2*pi*rand(1,num_Dian);
% �õ��������
x=x0+r.*cos(seta);
y=y0+r.*sin(seta);
x;
y;

% ������
 %plot(x,y,'*');
