function result = recurse_xor( bits )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

num_bits = length(bits);

if num_bits > 2
    result = recurse_xor(bits(1:num_bits-1));
    result = bitxor(result, bits(num_bits));
elseif num_bits == 1
    %error
    result = -1;
else
    result = bitxor(bits(1),bits(2));
end

