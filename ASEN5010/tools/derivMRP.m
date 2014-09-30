function sigma_dot = derivMRP( MRP, omega_body )
fcnPrintQueue(mfilename('fullpath'))
sigma_dot=0.25*((1-dot(MRP,MRP))*eye(3) + 2*vecSkew(MRP) + 2*(MRP*MRP'))*omega_body;

end

