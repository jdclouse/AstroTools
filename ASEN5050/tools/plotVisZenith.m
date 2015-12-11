function plotVisZenith(az, el)
%plotVisZenith Make a nice az/ze plot!
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% To go from Matlab's orientation to a clockwise position with north being
% up, need to reflect about 45 degree line (bisecting first and 3rd
% quadrants)

% Reflect about 45 degrees
fixed_az = pi/2 - az;

% Do the plot
polar(fixed_az, 90 - el * 180/pi)

% Grab the handle for the text. Set up the vectors of current, temp, and
% final labels
text = findall(gca, 'type', 'text');
text_angles = 0:30:330;
tmp_text = 'abcdefghijkl';
new_text_angles = [90:-30:0 330:-30:120];
% Set them to a temporary var to avoid redundantly and wrongly setting some
% labels
for ang = text_angles
    hText = findall(text, 'string', num2str(ang));
    set(hText, 'string', tmp_text(find(text_angles == ang,1)))
end
% Set them to the real thing
for ang = tmp_text
    hText = findall(text, 'string', num2str(ang));
    set(hText, 'string', num2str(new_text_angles(find(tmp_text == ang,1))))
end
xlabel('Azimuth (deg), Zenith (deg)')