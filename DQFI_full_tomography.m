function [bound_opt, V_opt, rho_opt]=DQFI_full_tomography(P0)
% Computes Detector SLD CR-bound for full POVM tomography.
%
% [c_sld, J, Lsld ] = sld_fun(S0,{S1,S2,S3})
% 
% P0 is the POVM with m-elements P0{1} through P0{m}
% Pn is the array (m x n) of POVM derivatives with Pn{j,k} the derivative
% of the j-th POVM element wrt the k-th parameter
% 
%


% requires YALMIP/Mosek


dim = size(P0{1},1);
mout = length(P0);

n = dim*dim;
m=mout;

% set derivatives for GMM tomography (independent of true values)
rhods = cell(m, 1);
rhods{m}=cell(1, n*(m-1));
for k=1:m-1
    rhods{k} = cell(1, n*(m-1));
    for k1=1:m-1
        for i=0:dim-1
            for j=0:dim-1
                if k1 == k
                    if i==0 && j==0
                        c = 1/sqrt(dim); % our normalisation for eye(dim) is c = 1 but DD use 1/sqrt(dim)
                    else
                        c = 1/sqrt(2);
                    end
                    rhods{k}{(k1-1)*n+i*dim+j+1} = c*GenGellMann(i,j,dim);
                    rhods{m}{(k1-1)*n+i*dim+j+1} = -c*GenGellMann(i,j,dim);
                else
                    rhods{k}{(k1-1)*n+i*dim+j+1} = zeros(dim);
                end
            end
        end
    end
end


% generate some random measurements
% POVMj = RandomPOVM(dim,m);
% 

[bound_opt, V_opt, rho_opt] = detCRB(P0, rhods);

disp(bound_opt)
%disp(VoptTomo)
disp(rho_opt)

end