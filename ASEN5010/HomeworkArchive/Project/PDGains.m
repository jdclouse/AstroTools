zeta=0.98;
Ts=70; %s
I11=14;
I22=16;
tolerance=0.02; % response at Ts will be within 2% of desired output
wn=-log(tolerance*sqrt(1-zeta*zeta))/zeta/Ts;

Kp1=wn*wn*I11;
Kd1=zeta*2*sqrt(I11*Kp1);
Kp2=wn*wn*I22;
Kd2=zeta*2*sqrt(I22*Kp2);

K=[Kp1 Kd1; Kp2 Kd2]