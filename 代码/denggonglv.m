clear all;

BS_Num = 1; 
BS_Radius = 100; 
UE_Num=20; %YONGHUSH
ZXD=10;
ZXD_UE=UE_Num/ZXD;
B= 3000000;%3MHz
B1=B/ZXD;
MaxBS_Power=35;%单位是dbm
MaxBS_Power1=10^((MaxBS_Power-30)/10);%单位是w
%Power=5:3:MaxBS_Power;
%Power1=10.^((Power-30)/10);%
max=8;
Power=0:1:max;
Iter_MonteCarlo_Max=100;
ThermalNoisedBm=-174; %-174dBm/Hz
noise=10^((ThermalNoisedBm-30)/10);%w 
pc=30;%单位是dbm
Pc=10^((pc-30)/10);%单位是w
Rmin=1;%单位是bit/s;
a=3;%损耗指数
QoS=Rmin;
% SubChannel_PowerRatio = zeros(ZXD,1);
M_array =[2,2,2,2,2,2,2,2,2,2];

EnergyEfficiency_NOMA= zeros(1,length(Power));
for Num_BS_Power=1:length(Power)
   
   % BSPower1=Power(Num_BS_Power);
    %BS_Power=10^((BSPower1-30)/10);
   BS_Power=Power(Num_BS_Power);
    Pn=BS_Power/ZXD;
     %  Gain1=sort(Gain,'ascend');
      for Iter_MonteCarlo = 1:Iter_MonteCarlo_Max
 BS.location_x = 0;  
        BS.location_y = 0;        
        [x,y]=DrawPoint( BS.location_x,BS.location_y,BS_Radius,UE_Num);
 for j=1:UE_Num
         UE(j).location_x=x(j);
        UE(j).location_y=y(j);
         UE(j).location_x;
         UE(j).location_x;         
end
        pathloss_BS_UE = zeros(1,UE_Num);
        distance_BS_UE = zeros(1,UE_Num);
for j=1:UE_Num
            distance_BS_UE(j) = compute_distance1(BS,UE(j));
            pathloss_BS_UE(j) = compute_pathloss_BS_UE_outdoor(a,distance_BS_UE(j));
end
       pathloss_BS_UE_SubChannel = pathloss_BS_UE;
        
BS_UE_Rayleigh0=1/sqrt(2)*randn(1,UE_Num)+1/sqrt(2)*randn(1,UE_Num)*1j;
        BS_UE_Rayleigh = BS_UE_Rayleigh0.*conj(BS_UE_Rayleigh0);
         Gain_BS_UE_SubChannel= zeros(1,UE_Num);
        for j=1:UE_Num
       Gain_BS_UE_SubChannel(j)=BS_UE_Rayleigh(j) /pathloss_BS_UE_SubChannel(j);
        end
       Gain1= Gain_BS_UE_SubChannel;
      fencu1=[Gain1(1),Gain1(2)];
 fencu2=[Gain1(3),Gain1(4)];
 fencu3=[Gain1(5),Gain1(6)];
 fencu4=[Gain1(7),Gain1(8)];
 fencu5=[Gain1(9),Gain1(10)];
 fencu1=[Gain1(11),Gain1(12)];
 fencu2=[Gain1(13),Gain1(14)];
 fencu3=[Gain1(15),Gain1(16)];
 fencu4=[Gain1(17),Gain1(18)];
 fencu5=[Gain1(19),Gain1(20)];
 h=[fencu1;fencu2;fencu3;fencu4;fencu5;fencu1;fencu2;fencu3;fencu4;fencu5;];
  h1=sort(h,2);
  arfa1=[];
  EE=[];
   for index=1:1:length(M_array)
    index = index
    M = M_array(index);
    Q = QoS*ones(1,M);
    Pn1=Pn*ones(1,M);

      nh1=[Gain1(2*index-1),Gain1(2*index)];
     nh=sort(nh1);
     
            ropt = zeros(1,M);
            A = (2.^(Q)-1)./(nh.*Pn1);
           for m = 1 : 1 : M-1
                ropt(m) = A(m)*(Pn*nh(m)*(1-sum(ropt(1:m-1))) + noise*B1)/(2^(QoS)); 
            end
            ropt(M) = 1-sum(ropt(1:M-1));
         
      EE(index)=[B1*log2( nh(2)*ropt(2)*Pn/(B1*noise))+B1*log2(nh(1)*Pn*ropt(1)/(Pn*ropt(2)*nh(1)+B1*noise))]/(Pn+Pc);
      
%arfa1;
 %  arfa=[rot90(arfa1)];
%arfa1=ropt(M);
%arfa=[rot90(ropt(M))];
 
 
   end
       EE;
   EnergyEfficiency_NOMA(Num_BS_Power)=EnergyEfficiency_NOMA(Num_BS_Power)+ sum(EE);
   %  AverageRs(index,:) = nanmean(Rs);
      end
end
    EnergyEfficiency_NOMA;
  EnergyEfficiency_NOMA=EnergyEfficiency_NOMA./(Iter_MonteCarlo_Max);
EnergyEfficiency_NOMA(1)=0;

%Power=5:3:MaxBS_Power;
plot(Power, EnergyEfficiency_NOMA,'r-o');
grid on;
legend('NOMA-DC')
ylabel('Energy effiency of system(bits/Joule)');
xlabel('Power of BS ');

