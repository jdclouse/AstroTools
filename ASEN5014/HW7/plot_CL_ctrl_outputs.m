function u = plot_CL_ctrl_outputs(r,F,K,y,t)

figure('Position',[0 0 1120 840])
u = [];
for ii = 1:length(y(:,1))
    u = [u; (F*r-K*y(ii,:)')'];
end
plot(t,u(:,1),t,u(:,2))
hold on
plot(t,pi/2*ones(1,length(t)),'r--')
plot(t,-pi/2*ones(1,length(t)),'r--')
legend('Aileron Angle', 'Rudder Angle')
xlabel('Time (s)')
ylabel('Deflection Angle (rad)')