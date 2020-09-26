function Q_p = EnergyEfficiencyfunc(Gain1,Gain2,beta,Pc,P,Power_k,noise,B1)
ThermalNoisedBm=-174; %-174dBm/Hz
noise=10^((ThermalNoisedBm-30)/10);%w 

B = beta; %  subchannel_num *1 matrix
Pk =Power_k;
P_n = Pc+Power_k;%  subchannel_num *1 matrix
F_p=zeros(length(B),1);
G_pk=zeros(length(B),1);
dG_pk=zeros(length(B),1);
for i=1:length(beta)
     if(Gain1(i)<Gain2(i))
         F_p(i)= (-log2(1+B(i).*P(i).*Gain2(i)/(noise*B1))-log2(noise*B1+P(i).*Gain1(i)))./(Pc+P(i));
         G_pk(i) =-log2(noise+B(i).*Pk(i).*Gain1(i))./(P_n(i)+Pc);
         %dG_pk(i) = ((log2(1+B(i).*Pk(i).*Gain1(i))-P_n(i).*B(i).*Gain1(i)./((1+B(i).*Pk(i).*Gain1(i)).*log(2)))./P_n(i).^2);
         dG_pk(i)=(log2(noise*B1+B(i).*Pk(i).*Gain1(i)))./(P_n(i)+Pc).^2-B(i).*Gain1(i)./((P_n(i)+Pc).*(noise*B1+B(i).*Pk(i).*Gain1(i)).*log(2));
     else
         F_p(i)= (-log2(1+B(i).*P(i).*Gain1(i))-log2(1+P(i).*Gain2(i)))./(Pc+P(i));
         G_pk(i) =-log2(1+B(i).*Pk(i).*Gain2(i))./(P_n(i));
         %dG_pk(i) = ((log2(1+B(i).*Pk(i).*Gain2(i))-P_n(i).*B(i).*Gain2(i)./((1+B(i).*Pk(i).*Gain2(i)).*log(2)))./P_n(i).^2);
         dG_pk(i)=(log2(noise+B(i).*Pk(i).*Gain2(i)))./(P_n(i)).^2-B(i).*Gain2(i)./(P_n(i).*(1+B(i).*Pk(i).*Gain2(i)).*log(2));
     end;

%dG_Pk =sum(log2(1+B.*Pk.*Gain2)-P_n.*B.*Gain2./((1+B.*Pk.*Gain2).*log2))./(P_n.^2);
%Q_p = Q_p+ F_p-G_Pk-dG_Pk.*(P(i)-Pk(i));
end   
Q_p = sum(F_p)-sum(G_pk)-sum(dG_pk'*(P-Pk));