clear

B = [sym('B0');sym('B1');sym('B2');sym('B3')];
Bp = [sym('Bp0');sym('Bp1');sym('Bp2');sym('Bp3')];
Bpp = [sym('Bpp0');sym('Bpp1');sym('Bpp2');sym('Bpp3')];

syms B0 B1 B2 B3 Bp0 Bp1 Bp2 Bp3 Bpp0 Bpp1 Bpp2 Bpp3

DCM_B = [B0^2 + B1^2 - B2^2 - B3^2, 2*(B1*B2 + B0*B3), 2*(B1*B3 - B0*B2);...
    2*(B1*B2 - B0*B3), B0^2 - B1^2 + B2^2 - B3^2, 2*(B2*B3 + B0*B1);...
    2*(B1*B3 + B0*B2), 2*(B2*B3 - B0*B1), B0^2 - B1^2 - B2^2 + B3^2]

DCM_Bp = [Bp0^2 + Bp1^2 - Bp2^2 - Bp3^2, 2*(Bp1*Bp2 + Bp0*Bp3), 2*(Bp1*Bp3 - Bp0*Bp2);...
    2*(Bp1*Bp2 - Bp0*Bp3), Bp0^2 - Bp1^2 + Bp2^2 - Bp3^2, 2*(Bp2*Bp3 + Bp0*Bp1);...
    2*(Bp1*Bp3 + Bp0*Bp2), 2*(Bp2*Bp3 - Bp0*Bp1), Bp0^2 - Bp1^2 - Bp2^2 + Bp3^2]
DCM_Bpp = [Bpp0^2 + Bpp1^2 - Bpp2^2 - Bpp3^2, 2*(Bpp1*Bpp2 + Bpp0*Bpp3), 2*(Bpp1*Bpp3 - Bpp0*Bpp2);...
    2*(Bpp1*Bpp2 - Bpp0*Bpp3), Bpp0^2 - Bpp1^2 + Bpp2^2 - Bpp3^2, 2*(Bpp2*Bpp3 + Bpp0*Bpp1);...
    2*(Bpp1*Bpp3 + Bpp0*Bpp2), 2*(Bpp2*Bpp3 - Bpp0*Bpp1), Bpp0^2 - Bpp1^2 - Bpp2^2 + Bpp3^2]

Composite_DCM = simplify(DCM_Bpp*DCM_Bp);
% solve(
