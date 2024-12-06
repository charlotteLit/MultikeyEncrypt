function [bin2_8] = decimal_binary(value,b)

bin2_8 = dec2bin(value)-'0';
if length(bin2_8) < b
    len = length(bin2_8);
    B = bin2_8;
    bin2_8 = zeros(1,b);
    for i=1:len
        bin2_8(b-len+i) = B(i); %不足8位前面补充0
    end 
end