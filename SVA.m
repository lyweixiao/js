M=4;
a=3;
dm=80;
sigma=-70;%dBm
pc=30;
h = 1/sqrt(2)*randn(1,M)+1/sqrt(2)*randn(1,M)*1j;
h = h.*conj(h);
nh = sort(h,'descend')./(1+dm.^a);  
l=0;
p1=4;%µ¥Î»ÊÇw;
p2=4;
p3=4;
p4=4;
r1=p1.*(h(1).^2)/(sigma);
r2=p2.*(h(4).^2)/(p1.*h(4).^2+sigma);
r3=p3.*(h(2).^2)/(sigma);
r4=p4.*(h(3).^2)/(p3.*h(3).^2+sigma);
a1=r1/(1+r1);
b1=log2(1+r1)-(r1/(1+r1)).*log2(r1);
a2=r2/(1+r2);
b2=log2(1+r2)-(r2/(1+r2)).*log2(r2);
a3=r3/(1+r3);
b3=log2(1+r3)-(r3/(1+r3)).*log2(r3);
a4=r4/(1+r4);
b4=log2(1+r4)-(r4/(1+r4)).*log2(r4);
for l=1:10;
    m=0;
    lamed=0;%
    dert=0.001;
    F=0.01;
    n=4;
    while F>dert;
        cvx_begin quiet
        variable q(n)
        maximize (a1*log(2.^q(1).*(h(1).^2)/log(2))+b1+ a2*log(2.^q(2).*(h(4).^2))/log(2)+b2+ a3*log(2.^q(3).*(h(2).^2))/log(2)+b3+a4*log(2.^q(4).*(h(3).^2))/log(2)+b4...
        -a1.*log(sigma)/log(2)-a2.*log(2.^q(1).*(h(4).^2)+sigma)/log(2) -a3.*log(sigma)/log(2)-a4.*log(2.^q(3).*(h(3).^2)+sigma))/log(2)-lamed*(p(1)+p(2)+p(3)+p(4)+pc)
        subject to
        log(1+2.^q(1).*(h(1).^2)/(sigma))/log(2)>=Rmin;
        log(1+2.^q(2).*(h(4).^2)/(2.^q(1).*h(4).^2+sigma))/log(2)>=Rmin;
        log(1+2.^q(3).*(h(2).^2)/(sigma))/log(2)>=Rmin;
        log(1+2.^q(4).*(h(3).^2)/(2.^q(3).*h(3).^2+sigma))/log(2)>=Rmin;
        q(1)+q(2)+q(3)+q(4)<=Pmax;
        q(1)>=0;
        q(2)>=0;
        q(3)>=0;
        q(4)>=0;
        cvx_end
     F=(a1*log2(2.^q(1).*(h(1).^2))+b1+ a2*log2(2.^q(2).*(h(4).^2))+b2+ a3*log2(2.^q(3).*(h(2).^2))+b3+a4*log2(2.^q(4).*(h(3).^2))+b4...
        -a1.*log2(sigma)-a2.*log2(2.^q(1).*(h(4).^2)+sigma) -a3.*log2(sigma)-a4.*log2(2.^q(3).*(h(3).^2)+sigma))-lamed*(p(1)+p(2)+p(3)+p(4)+pc);
    lamed=(a1*log2(2.^q(1).*(h(1).^2))+b1+ a2*log2(2.^q(2).*(h(4).^2))+b2+ a3*log2(2.^q(3).*(h(2).^2))+b3+a4*log2(2.^q(4).*(h(3).^2))+b4...
        -a1.*log2(sigma)-a2.*log2(2.^q(1).*(h(4).^2)+sigma) -a3.*log2(sigma)-a4.*log2(2.^q(3).*(h(3).^2)+sigma))/(p(1)+p(2)+p(3)+p(4)+pc);
    m=m+1;
    end
    p1=2.^q(1);
    p2=2.^q(2);
    p3=2.^q(3);
    p4=2.^q(4);
    r1=p1.*(h(1).^2)/(sigma);
    r2=p2.*(h(4).^2)/(p1.*h(4).^2+sigma);
    r3=p3.*(h(2).^2)/(sigma);
    r4=p4.*(h(3).^2)/(p3.*h(3).^2+sigma);
    a1=r1/(1+r1);
    b1=log2(1+r1)-(r1/(1+r1)).*log2(r1);
    a2=r2/(1+r2);
    b2=log2(1+r2)-(r2/(1+r2)).*log2(r2);
    a3=r3/(1+r3);
    b3=log2(1+r3)-(r3/(1+r3)).*log2(r3);
    a4=r4/(1+r4);
    b4=log2(1+r4)-(r4/(1+r4)).*log2(r4);
    l=l+1;
end




