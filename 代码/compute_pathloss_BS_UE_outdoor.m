function  [pathloss] = compute_pathloss_MBS_MUE_outdoor(b,mDistance)
%pathloss_dB=sqrt(1+mDistance^0.1);
pathloss=1+mDistance^b;

%pathloss_dB = 133.6+35*log10(mDistance/1000);%pathloss dB model
%MUEAntennaGain = 0; 
%MUEShadowFading=8; %eNB-UE log normal shadowing 10
%UETotalPathLossIndB = MUEAntennaGain-MUEShadowFading-pathloss_dB;
%pathloss=10^(UETotalPathLossIndB/10);%Total path