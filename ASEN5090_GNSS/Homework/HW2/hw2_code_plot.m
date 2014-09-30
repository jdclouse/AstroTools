function hw2_code_plot(code)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

figure
subplot(2,1,1)
stairs(code.CA_code(1:17)), ylim([-0.5, 1.5])
title(sprintf('PRN %i: First 16 bits', code.PRN_Number))
xlabel('Bit')
ylabel('Value')
subplot(2,1,2)
stairs(code.CA_code(end-16:end)), ylim([-0.5, 1.5])
title(sprintf('PRN %i: Last 16 bits', code.PRN_Number))
xlabel('Bit')
ylabel('Value')