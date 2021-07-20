function [theta] = theta_rotate_1(lenu,lenv,max_x,max_y,rotate)
theta = zeros(lenu,lenv);

for i = 1:lenu
    for j = 1:lenv
        x=i-max_y;
        y=j-max_x;
%         old_ori = rot_map(i,j);
        r = sqrt(x^2+y^2);
        if (r<10) % size of guiding sector
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
%             rot_map(i,j) = rotate;
        end
    end
end

for i = max_y-1:max_y+1
    for j = max_x-1:max_x+1
%         x_dist = i-max_x;
%         y_dist = j-max_y;
%         r = sqrt(x_dist^2+y_dist^2);
%         if (r<=1)
            theta(i,j) = 3;
%         end
    end
end

theta(abs(theta)<=2.8) = 0;
theta(abs(theta)>2.8) = 1;
