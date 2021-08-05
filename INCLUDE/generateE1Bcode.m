function E1Bcode = generateE1Bcode(PRN, BOC)

codes = load('codes_E1B.mat');
E1B=codes.codes_E1B(:,PRN);

chips=4092;
if BOC==1
    E1Bcode=ones(1,4*chips);
    sub=ones(1,4*chips);
    
    for i=1:(4*chips)
        if rem(i,2) == 0
            sub(i)= -1;
        end %if
            E1Bdouble(i) = E1B(ceil(i/4));
            E1Bcode(i) = E1Bdouble(i)*sub(i);
    end %for
else 
   E1Bcode = E1B; 
end %if

    