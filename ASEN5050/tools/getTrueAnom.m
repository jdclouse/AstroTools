function f = getTrueAnom( r,v,mu,t)
%getTrunAnom get only trun anom from r,v,mu
%   Output in radians
fcnPrintQueue(mfilename('fullpath'))
[~,~,~,~,~,f] = cart2OE(r,v,mu);

global progress;
if f*180/pi > progress(1)
    fprintf('Progress: %d degrees\n',progress(1))
    fprintf('Day %.1f\n',t/3600/24)
    if length(progress) > 1
    progress = progress(progress > f*180/pi);
    end
    
end

end

