function [b, temp_mat] = lip_param_2(mu,x,discrep, exp_size, bin_param)
% mu is the vector of mu_hat, is supposed to have dimension k * 1
% x is a matrix of systems, is supposed to have dimension k * n
% discrep is the vector of discrepency, which has dimension k * 1
% bin_param is the parameter for binary search
    % exp_size = 29966;
    a = 0;
    b = 1000;
    while((b-a) > bin_param)
        flag = 0;
        t = (a+b)/2;
        for i=1:exp_size
            mat_mu = pdist2(mu(i), mu);
            mat_dis = discrep(i) + discrep;
            mat_x = pdist2(x(i,:), x, 'chebychev');
            diff_mat = mat_dis - mat_mu;
            temp_mat = t*mat_x + diff_mat;
            if min(temp_mat,[],"all")<0
                a = t;
                flag = 1;
                break
            end
        end
        if flag==0
            b = t;
        end
        disp(b)
        disp(a)
    end
end
