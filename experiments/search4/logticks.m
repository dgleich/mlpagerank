function xs = logticks(x)
% LOGTICKS Return logspace marks at powers of 10 up to x
xseq = [1 ];
xs = [];
power = 1;
while 1
    for i=1:numel(xseq)
        if xseq(i)*power > x
            break;
        end
        xs(end+1) = xseq(i)*power;
    end
    if xseq(i)*power > x
        break;
    end
    power = power*10;
end