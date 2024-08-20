function b = lip_param(mu,x,discrep, bin_param)
% mu is the vector of mu_hat, is supposed to have dimension k * 1
% x is a matrix of systems, is supposed to have dimension k * n
% discrep is the vector of discrepency, which has dimension k * 1
% bin_param is the parameter for binary search
    mat_mu = pdist2(mu, mu);
    mat_dis = discrep + discrep';
    mat_x = pdist2(x, x, 'chebychev');
    a = 0;
    b = 100000;

    diff_mat = mat_dis - mat_mu;

    while(b-a >bin_param)
        t = (a+b)/2;
        temp_mat = t*mat_x + diff_mat;
        temp_mat = (temp_mat<0);
        if sum(temp_mat,"all")>0
            a = t;
        else
            b = 5;
        end
    end
end