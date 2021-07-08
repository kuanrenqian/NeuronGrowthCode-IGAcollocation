function [theta] = theta_rotate(theta,lenu,lenv,max_x,max_y,rotate)
for i = 1:lenu
    for j = 1:lenv
        x=i-max_y;
        y=j-max_x;
        theta(i,j) = (atan2(y,x));
        if isnan(theta(i,j))
            theta(i,j) = 1;
        end
        theta(i,j) = theta(i,j)+rotate;
        if(theta(i,j)>pi)
            theta(i,j) = theta(i,j) - 2*pi;
        elseif (theta(i,j)<-pi)
            theta(i,j) = theta(i,j) + 2*pi;
        end
    end
end
val = 3;
for i = max_y-1:max_y+2
    for j = max_x-1:max_x+1
        theta(i,j) = val;
        val = -val;
    end
end
