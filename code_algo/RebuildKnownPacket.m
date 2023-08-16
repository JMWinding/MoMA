function y = RebuildKnownPacket(y, dchip, hp, delta)
ylen = length(y);

yd = conv(hp,dchip);
ydlen = length(yd);

d1 = max(0,delta);
d2 = max(0,ylen-ydlen-delta);
y(1+d1:ylen-d2) = y(1+d1:ylen-d2) ...
    + yd(1-delta+d1:ylen-delta-d2);
    
end