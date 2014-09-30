function out = normalized_correlation( in_code1, in_code2, shift )
%normalized_correlation Output the correlation between two codes.
%codes are row vectors
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

%TODO Error check the codes
code_length = length(in_code1);

%Set the bits from 0,1 to 1,-1
code1 = in_code1*-1;
code2 = in_code2*-1;

for ii = 1:code_length
    if code1(ii) == 0
        code1(ii) = 1;
    end
    if code2(ii) == 0
        code2(ii) = 1;
    end
end

%Shift the second code
if shift > 0
    code2 = [code2(shift+1:end), code2(1:shift)];
elseif shift < 0
    code2 = [code2(code_length+shift+1:code_length), code2(1:code_length+shift)];
end

%Sum and divide
total = 0;
for ii = 1:code_length
    total = total + code1(ii) * code2(ii);
end

out = total/code_length;

end

