% moving x inside (lb,ub) interval
function xx=move_inside_period(x,lb,ub)

xx=x;
[nrow,ncol]=size(xx);
for i=1:nrow
    if (xx(i)<lb(i))|(xx(i)>ub(i))
        ord=floor((xx(i)-lb(i))/(ub(i)-lb(i)));
        ord=mod(ord,2); % f.e. mod(13,5) = 3
        xx(i)=mod(xx(i)-lb(i),ub(i)-lb(i))+lb(i);
        if (xx(i)<lb(i))
            xx(i)=xx(i)+ub(i)-lb(i);
        end
        if ord==1
            xx(i)=lb(i)+ub(i)-xx(i);
        end
    end
end
