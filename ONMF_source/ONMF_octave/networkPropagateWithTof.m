function new_mat = networkPropagateWithTof(src_mat, adj_mat, alpha, max_iter)

    %alpha = .6;
    itcnt = 0;
    tof = 1e-4;
%   max_iter = 1500;

    prev_src_mat = src_mat; 
    while (itcnt < max_iter)        
        new_mat = alpha * prev_src_mat * adj_mat + (1-alpha) * src_mat;        
       
        l1residual = norm(prev_src_mat - new_mat, 1);
        prev_src_mat = new_mat;
        
        if (l1residual == 0)
            fprintf(1,'Converged %d\n', itcnt);
            
            break;
        end
        
        itcnt = itcnt + 1;
    end
    
    fprintf(1, '%d-th iter: l1residua\nl', itcnt);
