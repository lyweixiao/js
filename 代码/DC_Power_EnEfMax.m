function [PowerAllocation_SubChannel,  EnergyEfficiency] = DC_Power_EnEfMax(Pn,Pc,MaxBS_Power1,ZXD, SubChannel_MatchGain1,SubChannel_MatchGain2, arfa)

beta = arfa; %  subchannel_num *1 matrix
Gain1 = SubChannel_MatchGain2;%  subchannel_num *1 matrix
Gain2 = SubChannel_MatchGain1;%  subchannel_num *1 matrix

   % P_initial = SubChannel_SumRateRatio.*BS_Power;
  P_initial = Pn.*ones(ZXD,1);
%  P_initial = zeros(SubChannel_Num,1);
   A=-eye(ZXD);
   b=zeros(ZXD,1);
   Aeq=ones(1,ZXD);
   beq =MaxBS_Power1;
   Power_k = P_initial;
%   P = sym('P',[SubChannel_Num,1]);
   % lb = (BS_Power/(SubChannel_Num*2)).*ones(1,SubChannel_Num)
  lb = (MaxBS_Power1/(3/2*ZXD)).*ones(1,ZXD);
 %  lb =0.*ones(1,ZXD);
   ub = MaxBS_Power1.*ones(1,ZXD);
Epsilon = 0.01;
%syms P;
Q2=10;
Q1=0;
k=0;
while(k<10 )
   if(abs(Q2-Q1)<Epsilon)
      break;
   end;
    if(k>0)
        Q1=Q2;
    end
   % Qp=SumRate_DC(Gain1,Gain2,beta,P,Power_k)
    options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
    options.Algorithm = 'sqp';
   [P,Q2] = fmincon(@(P) EnergyEfficiencyfunc(Gain1,Gain2,beta,Pc,P,Power_k),P_initial,A,b,Aeq,beq,lb,ub,[],options);
   Power_k=P; 
   k=k+1
end
PowerAllocation_SubChannel = P;
 EnergyEfficiency = -Q2;
%EnergyEfficiency = -Q2./(P+Pc);