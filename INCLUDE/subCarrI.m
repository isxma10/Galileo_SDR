function subI = subCarrI(points)

subI=ones(1,2*points);
for i=1:(2*points)
    if rem(i,2) == 0
        subI(i)= -1;
    end %if
end %for


