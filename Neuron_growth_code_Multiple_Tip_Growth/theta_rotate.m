function [theta] = theta_rotate(lenu,lenv,Max_x,Max_y,initial_angle,rotate,size_Max)
theta = zeros(lenu,lenv);
ttt = zeros(lenu,lenv);
for l=1:size_Max
    max_x = Max_x(l);
    max_y = Max_y(l);
%     old_ori = initial_angle(l)+rotate(l);
    old_ori = initial_angle(l);
    if(old_ori>0)
        old_ori = pi-old_ori;
    else
        old_ori = -(pi+old_ori);
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
                ttt(i,j) = ttt(i,j)+old_ori;
                if(ttt(i,j)>pi)
                    ttt(i,j) = ttt(i,j) - 2*pi;
                elseif (ttt(i,j)<-pi)
                    ttt(i,j) = ttt(i,j) + 2*pi;
                end
                ttt(abs(ttt)<2.5) = 0;
            end
        end
    end
%     val = 3;
%     for i = max_y-5:max_y+5
%         for j = max_x-5:max_x+5
%             x_dist = i-max_x;
%             y_dist = j-max_y;
%             r = sqrt(x_dist^2+y_dist^2);
%             if (r<2)
%                 ttt(i,j) = val;
%                 val = -val;
%             end
%         end
%     end
    for i = max_y-1:max_y+1
        for j = max_x-1:max_x+1
            ttt(i,j) = 1;
%             val = -val;
        end
    end
    
    theta = theta+ttt;    
end


theta(abs(theta)>0) = 1;