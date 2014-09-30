function vectorPlot(time, vector_array)
[time_r, time_c] = size(time);
if time_r ~= 1
   fprintf('Time array num rows != 1. Returning.\n')
   return
end
[array_r, array_c] = size(vector_array);
if array_r ~= 3
   fprintf('bad\n')
   return
elseif array_c ~= time_c
  fprintf('Time array and vector array have differing sizes. Returning.\n')
  return
end
  
graphics_toolkit('gnuplot')
figure
hold on 
plot(time, vector_array(1,:))
plot(time, vector_array(2,:),'r')
plot(time, vector_array(3,:), 'g')