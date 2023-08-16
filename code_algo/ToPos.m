function rval = ToPos(input)
rval = input;
rval(rval < 0) = 0;
end