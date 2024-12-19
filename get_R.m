function [value, i] = get_R(input_matrix, i)
    len = length(input_matrix);
    if i == 0
        value = input_matrix(1,1);
    elseif i > 0
        value = input_matrix(1,i+1);
    else
        value = input_matrix(len, len+i);
    end
end