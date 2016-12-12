% Problem 1 part (i)
function [ret_int] = randint( min, max)

    ret_int = round((min - 0.5) + rand() * ((max+0.5) - (min-0.5)));

end

