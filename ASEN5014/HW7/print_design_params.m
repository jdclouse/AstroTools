function print_design_params(y,t,PO,PS,Ts,ii)

overshoot_under_lim = 'False';
if max(abs(y(:,ii)-.25)) < (1+PO)*.25
    overshoot_under_lim = 'True';
end
fprintf(['Max overshoot within ' num2str(PO*100) '%%: ' ...
    overshoot_under_lim ' (%.3f)\n'], max(y(:,ii)))
settled_under_lim = 'False';
if max(abs(y(t>=Ts,ii) - 0.25)) < (1+PS)*.25
    settled_under_lim = 'True';
end
fprintf(['Settled within ' num2str(Ts) ' sec: ' settled_under_lim '\n\n'])