function [bound_opt, Lsld, V_opt, rho_opt, trQFI, QFImat]=DQFI_multi_para(P0,Pn)
% Computes Detector SLD CR-bound for parameter estimation.
%
% [c_sld, J, Lsld ] = sld_fun(S0,{S1,S2,S3})
% 
% P0 is the POVM with m-elements P0{1} through P0{m}
% Pn is the array (m x n) of POVM derivatives with Pn{j,k} the derivative
% of the j-th POVM element wrt the k-th parameter
% 

% requires YALMIP/Mosek


dim = size(P0{1},1);
mout = length(P0);
n = length(Pn{1});

Lsld = cell(mout, n);

for i=1:mout
    [vec,lam]=eig(P0{i});
    lam=diag(lam);

    % SLD operator Lsld, intialise to zero
    for l=1:n
        Lsld{i,l} = zeros(dim);
    end

    % update SLD operator L
    for j=1:dim
        for k=1:dim
            if lam(j)+lam(k)==0
                continue;
            else
                vj=vec(:,j);
                vk=vec(:,k);
                for l=1:n
                    Lsld{i, l} = Lsld{i, l} + 2*vj*vj'*Pn{i}{l}*vk*vk'/(lam(j)+lam(k));
                end
            end
        end
    end
end

% Define optimization variables using YALMIP
V = sdpvar(n, n, 'symmetric', 'real');    % Real symmetric matrix V
rho = sdpvar(dim, dim, 'hermitian', 'complex');  % Complex Hermitian matrix rho
Q = sdpvar(n, n, 'symmetric', 'real');

% calculate Tr DQFI CRB for reference:
% Define Q_theta(rho) according to Theorem 5
Qtr = zeros(n, n);
for j = 1:n
    for k = 1:n
        % Initialize element (j,k) of Q
        Qtr_temp = 0;
        % Sum over all l in [m]
        for l = 1:mout
            % Calculate the two terms inside the trace
            term1 = Lsld{l,j} * P0{l} * Lsld{l,k};
            term2 = Lsld{l,k} * P0{l} * Lsld{l,j};
            % Add to Q_jk
            Qtr_temp = Qtr_temp + 0.5 * trace((term1 + term2));
        end
        Qtr(j,k) = Qtr_temp;
    end
end
trQFI = trace(inv(Qtr));

QFImat = Qtr;

% Define Q_theta(rho) according to Theorem 5
for j = 1:n
    for k = 1:n
        % Initialize element (j,k) of Q
        Q_jk = 0;
        % Sum over all l in [m]
        for l = 1:mout
            % Calculate the two terms inside the trace
            term1 = Lsld{l,j} * P0{l} * Lsld{l,k};
            term2 = Lsld{l,k} * P0{l} * Lsld{l,j};
            % Add to Q_jk
            Q_jk = Q_jk + 0.5 * trace((term1 + term2) * rho);
        end
        Q(j,k) = Q_jk;
    end
end

bigmat = [V, eye(n); eye(n), Q];


% Constraint: V ≽ (Q_theta(rho))^(-1)
% Using Schur complement to avoid explicit inversion:
% [V, I; I, Q] ≽ 0 is equivalent to V ≽ Q^(-1)
schur_constraint = [bigmat>=0];

% Additional constraints: rho ≽ 0, trace(rho) = 1
rho_constraint1 = [rho>=0];
rho_constraint2 = [trace(rho)==1];

% Combine all constraints
constraints = [schur_constraint,rho_constraint1,rho_constraint2];

% Objective: minimize trace(W*V)
objective = trace(V);

% Set options for MOSEK
options = sdpsettings('verbose',0,'debug',0,'solver','mosek');

% Solve the SDP
result = optimize(constraints, objective, options);

% Extract optimal values
if result.problem == 0
    V_opt = value(V);
    rho_opt = value(rho);
    bound_opt = value(objective);
else
    warning('Optimization failed with status: %s', result.info);
    V_opt = [];
    rho_opt = [];
    bound_opt = -1;
end