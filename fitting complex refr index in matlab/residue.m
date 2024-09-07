% function to calculate residue, is used in fminsearch
function f=residue(x,Data,terms_num,terms_sz,arg,lb,ub,inum)

f=0;

[Data_nrow,Data_ncol]=size(Data);
for i=1:Data_nrow
    % epsilon calculated with current fitting
    eps=calculate_eps(x,Data(i,1),terms_num,terms_sz,arg,lb,ub);
    % difference between fitting and tabular data
    fre=real(eps)-Data(i,2);
    fim=imag(eps)-Data(i,3);

%    fim=10*log(abs(imag(eps)-Data(i,3)))/log(10);     
%    fre=fre/Data(i,2);
%    fim=fim/Data(i,3);

%    use this block to get accurate fitting for small imaginary epsilon values
    mult=2;
    if Data(i,3)<0.05         
       fim=fim*100*mult;
    elseif Data(i,3)<0.5
       fim=fim*10*mult;
    elseif Data(i,3)<5
       fim=fim*1*mult;
    else
       fim=fim*1*mult;
    end

    % calculating residue
    f=f+fre*fre+fim*fim;
end
f=sqrt(f);
