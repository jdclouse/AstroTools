function abs_exponent = dec_places( input )
%dec_places Return the sci-no exponent of numbers less than 1. It's quick
%and dirty.
fcnPrintQueue(mfilename('fullpath')); % Add this code to code app 

abs_exponent = 0;
comp_value = 1;
while input < comp_value
    abs_exponent = abs_exponent+1;
    comp_value = comp_value/10;
end
abs_exponent;
end

