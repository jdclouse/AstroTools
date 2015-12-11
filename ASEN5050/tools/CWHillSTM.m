function STM = CWHillSTM(w,t)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

s = sin(w*t);
c = cos(w*t);
x_t = [4-3*c, 0, 0, s/w, -2/w*c+2/w, 0];
y_t = [6*s-6*w*t, 1, 0, 2/w*c-2/w, 4/w*s-3*t, 0];
z_t = [0, 0, c, 0, 0, s/w];
xd_t = [3*w*s, 0, 0, c, 2*s, 0];
yd_t = [6*w*c-6*w, 0, 0, -2*s, 4*c-3, 0];
zd_t = [0, 0, -w*s, 0, 0, c];

STM = [x_t; y_t; z_t; xd_t; yd_t; zd_t];