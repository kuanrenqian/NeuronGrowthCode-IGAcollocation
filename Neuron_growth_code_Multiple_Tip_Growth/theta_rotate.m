function [theta] = theta_rotate(lenu,lenv,Max_x,Max_y,angle,size_Max)
theta = zeros(lenu,lenv);
ttt = zeros(lenu,lenv);
for l=1:size_Max
    max_x = Max_x(l);
    max_y = Max_y(l);
    if(angle(l)>0)
        angle(l) = pi-angle(l);
    else
        angle(l) = -(pi+angle(l));
    end
    for i = 1:lenu
        for j = 1:lenv
            y=i-max_y;
            x=j-max_x;
            r = sqrt(x^2+y^2);
            if (r<10) % size of guiding sector
                ttt(i,j) = atan2(x,y);
                if isnan(ttt(i,j))
                    ttt(i,j) = 1;
                end
                ttt(i,j) = ttt(i,j)+angle(l);
                if(ttt(i,j)>pi)
                    ttt(i,j) = ttt(i,j) - 2*pi;
                elseif (ttt(i,j)<-pi)
                    ttt(i,j) = ttt(i,j) + 2*pi;
                end
                ttt(abs(ttt)<2.5) = 0;
            end
        end
    end
    
    for i = max_y-1:max_y+1
        for j = max_x-1:max_x+1
            ttt(i,j) = 1;
        end
    end
    theta = theta+ttt;    
end
theta(abs(theta)>0) = 1;