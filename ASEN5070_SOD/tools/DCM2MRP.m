function MRP = DCM2MRP( DCM ) 
%DCM2MRP Turn DCM to MRP vector
fcnPrintQueue(mfilename('fullpath'))
zeta=sqrt(trace(DCM)+1);
sig_skew = (DCM'-DCM)/(zeta*(zeta+2));
MRP = [sig_skew(3,2);sig_skew(1,3);sig_skew(2,1)];
end

