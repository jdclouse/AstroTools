function angles_out = within180(angle_set)
fcnPrintQueue(mfilename('fullpath'))
for i = 1:3
   if angle_set(i) > 180
      angle_set(i) -= 360;
   elseif angle_set(i) < -180
      angle_set(i) += 360;
   end
end
angles_out = angle_set;
end