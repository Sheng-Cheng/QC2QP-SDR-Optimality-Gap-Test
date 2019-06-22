function status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,epsilon2,verbosity)
% The QC2QP_SDR_optimalityGap_test is a function that can test whether a 
% type of optimization problem, named quadratic program with two quadratic
% constraints (QC2QP), can be solved by its semidefinite relaxation, in 
% which case we say there is no optimality gap. Specifically, taking an 
% example of a two dimensional QC2QP, the problem has the following form:
%
% minimize    z^T*[a0 b0 b0 c0]*z + 2*[d0 e0]^T*z + f0
% {z \in R}
% subject to  z^T*[a1 b1 b1 c1]*z + 2*[d1 e1]^T*z + f1 <= 0
%             z^T*[a2 b2 b2 c2]*z + 2*[d2 e2]^T*z + f2 <= 0

% Please make sure CVX and YALMIP are both installed on your machine. They
% are required by this test.

% Input:
%   - M0: homogeneous form of paramters for the objective function. An
%  	      exemple for M0 in above problem is the following:
%         M0 = [f0 d0 e0;
%               d0 a0 b0;
%               e0 b0 c0];
%   - M1: homogeneous form of paramters for the first constraint function.
%         An exemple for M1 in above problem is the following:
%         M1 = [f1 d1 e1;
%               d1 a1 b1;
%               e1 b1 c1];
%   - M2: homogeneous form of paramters for the first constraint function.
%         An exemple for M1 in above problem is the following:
%         M2 = [f2 d2 e2;
%               d2 a2 b2;
%               e2 b2 c2];
%   - epsilon2: tolerance value for purification. This value is determined
%                by the user, while the default value is 1e-5.
%   - verbosity: switch for the amount of messages displayed. 0 for no
%                messages and 1 for the messages comming from CVX and this
%                test. The default value is 0.

% Output: a struct named 'status' which contains the following fields:
%   - epsilon2: to record the purification tolerance
%   - SlaterSP: to record whether Slater's condition holds for the
%               semidefinite relaxation (SP). 1 for holding and 0 for not.
%   - SlaterSD: to record whether Slater's condition holds for the
%               dual problem (SD). 1 for holding and 0 for not.
%   - noOptimalityGap: indicator of whether there exists an optimality gap.
%                      1 for no optimality gap, 0 for optimality gap, and 
%                      NaN for undetermined.
%   - exception.code: indicator of whether an exception appears.
%                     0 for no exception;
%                     1 for problem ill-posed;
%                     2 for rank(X_hat) >= 3;
%                     3 for error showing up when rank(Z_hat) < n-1 while
%                     rank(X_hat) == 2
%                     4 for alpha1 undetermined
%                     (Please email the author (cheng@terpmail.umd.edu)
%                     if the exception code is not 0)
%   - exception.M0: value of M0 that causes an
%                   exception. It only shows up when exception.code is not 0.
%   - exception.M1: value of M1 that causes an
%                   exception. It only shows up when exception.code is not 0.
%   - exception.M2: value of M2 that causes an
%                   exception. It only shows up when exception.code is not 0.
%   If the problem data yield Slater's condition holding for both (SP) and
%   (SD), then the following fields are attached.
%   - SP_optval: optimal value of (SP)
%   - SP_optsol: optimal solution of (SP), X_hat
%   - SP_x1: rank-one decomposition of X_hat
%   - SP_x2: rank-one decomposition of X_hat (this field only exists if rank(X_hat) == 2)
%   - SD_optval: optimal value of (SD)
%   - SD_optsol.Z: optimal solution of (SD), Z_hat
%   - SD_optsol.y0: optimal solution of (SD), y_hat0
%   - SD_optsol.y1: optimal solution of (SD), y_hat1
%   - SD_optsol.y2: optimal solution of (SD), y_hat2
%   - rank_X: rank of the solution to (SP) with tolerance epsilon2
%   - rank_Z: rank of the solution to (SD) with tolerance epsilon2
%   If the problem data yield no optimality gap, then the following fields are
%   attacehd.
%   - primal_optval: optimal value of the primal problem
%   - primal_optsol: optimal solution of the primal problem
%
% Sheng Cheng (cheng@terpmail.umd.edu), Nov. 2018.

if nargin <= 2
    error('Insufficient inputs.\n');
elseif nargin == 3
    epsilon2 = 1e-5; % default value of epsilon2 is 1e-5
    verbosity = 0; % default value of verbosity is 0
elseif nargin == 4
    verbosity = 0;
    if isempty(epsilon2) % default value of epsilon2 is 1e-5
        epsilon2 = 1e-5;
    end
elseif nargin == 5
    if (verbosity~=0) && (verbosity ~=1)
        error('Verbosity must be 0 or 1.\n');
    else
        if isempty(epsilon2)
            epsilon2 = 1e-5; % default value of epsilon2 is 1e-5
        end
    end
end

if ~issymmetric(M0)
    error('M0 is not symmetric. Abort.');
end

if ~issymmetric(M1)
    error('M1 is not symmetric. Abort.');
end

if ~issymmetric(M2)
    error('M2 is not symmetric. Abort.');
end

if ~isreal(M0)
    error('M0 is not real. Abort.');
end

if ~isreal(M1)
    error('M1 is not real. Abort.');
end

if ~isreal(M2)
    error('M2 is not real. Abort.');
end

% set the tag of using relative rank
useRelRankTolTag = 0;

% determine the dimension of problem
n = size(M0,1)-1;

if ~((size(M1,1) == n+1) && (size(M2,1) == n+1))
    error('M0, M1, and M2 are not of the same dimension. Abort.');
end

% initilize the output struct
status = struct;

status.epsilon2 = epsilon2;
status.SlaterSP = []; % 0 for failing and 1 for passing the Slater's condition check for (SP)
status.SlaterSD = []; % 0 for failing and 1 for passing the Slater's condition check for (SD)
status.noOptimalityGap = NaN; % NaN: existence of the optimality gap cannot be determined by this test because of assumption violation
% 0: there is no optimality gap
% 1: there is an optimality gap

I00 = diag([1 zeros(1,n)]);

SlaterCheck = 1; % 0 for failing the Slater's condition check, 1 for passing the Slater's condition check

noOptimalityGapTag = []; % 0 for optimality gap, nonzero for no optimality gap

% SlaterCheck: Slater's condition for (SP)
cvx_begin SDP quiet
cvx_solver SDPT3
variable X_hat(n+1,n+1) symmetric
subject to
trace(M1*X_hat) <= -epsilon2;
trace(M2*X_hat) <= -epsilon2;
trace(I00*X_hat) == 1;
X_hat >= epsilon2*eye(n+1);
cvx_end

% If the above problem is infeasible, then Slater's condition does not hold
% for (SP).
if cvx_optval == Inf
    fprintf('Slater condition does not hold for (SP).\nAbort.\n');
    SlaterCheck = SlaterCheck && 0;
    status.SlaterSP = 0;
else
    status.SlaterSP = 1;
end

% SlaterCheck: Slater's condition for (SD)
cvx_begin SDP quiet
cvx_solver SDPT3
variables y_hat0 y_hat1 y_hat2
subject to
M0-y_hat0*I00+y_hat1*M1+y_hat2*M2 >= epsilon2*eye(n+1);
y_hat1 >= epsilon2;
y_hat2 >= epsilon2;
cvx_end

% If the above problem is infeasible, then Slater's condition does not hold for (SD)
if cvx_optval == Inf
    fprintf('Slater condition does not hold for (SD).\nAbort.\n');
    SlaterCheck = SlaterCheck && 0;
    status.SlaterSD = 0;
else
    status.SlaterSD = 1;
end

if SlaterCheck
    
    status.SP_optval = []; % optimal value of (SP)
    status.SP_optsol = []; % optimal solution of (SP)
    
    status.SD_optval = []; % optimal value of (SD)
    status.SD_optsol = struct; % optimal solution of (SD)
    
    % solving (SP) and it dual (SD)
    switch verbosity
        case 1
            cvx_begin SDP
        case 0
            cvx_begin SDP quiet
    end
    cvx_solver SDPT3
    variable X_hat(n+1,n+1) symmetric
    dual variables y_hat0 y_hat1 y_hat2 Z_hat
    minimize (trace(M0*X_hat))
    subject to
    y_hat1:trace(M1*X_hat) <= 0;
    y_hat2:trace(M2*X_hat) <= 0;
    y_hat0:trace(I00*X_hat) == 1;
    Z_hat:X_hat >= 0;
    cvx_end
    
    if cvx_optval == -Inf
        % Since Slater's condition holds for (SD), the solution of (SP)
        % shouldn't be -Inf. We will switch to the backup solver, sedumi
        fprintf('SDPT3 cannot solve (SP), switching to Sedumi.\n');
        switch verbosity
            case 1
                cvx_begin SDP
            case 0
                cvx_begin SDP quiet
        end
        cvx_solver sedumi
        variable X_hat(n+1,n+1) symmetric
        dual variables y_hat0 y_hat1 y_hat2 Z_hat
        minimize (trace(M0*X_hat))
        subject to
        y_hat1:trace(M1*X_hat) <= 0;
        y_hat2:trace(M2*X_hat) <= 0;
        y_hat0:trace(I00*X_hat) == 1;
        Z_hat:X_hat >= 0;
        cvx_end
        
        if cvx_optval == -Inf
            fprintf('Problem ill-posed. \n')
            status.exception.code = 2; % 2 for problem ill-posed
            status.exception.M0 = M0;
            status.exception.M1 = M1;
            status.exception.M2 = M2;
            status.exception.epsilon2 = epsilon2;
            return % return the control to the invoking function
        end
    end
    
    status.SP_optval = cvx_optval;
    status.SP_optsol = X_hat;
    
    status.SD_optval = cvx_optval;
    status.SD_optsol.Z = Z_hat;
    status.SD_optsol.y0 = y_hat0;
    status.SD_optsol.y1 = y_hat1;
    status.SD_optsol.y2 = y_hat2;
    
    % eigen decompose Z_hat and X_hat
    [V_X,D_X] = eig(X_hat);
    [V_Z,D_Z] = eig(Z_hat);
    
    % compute the numerical rank of X with tolerance epsilon2
    rank_X = sum(abs(diag(D_X)) >= epsilon2);
    if rank_X >= 3
        % if absolute tolerance gives incorrect rank, then switch to
        % relative tolerance automatically
        fprintf('Absolute rank tolerance causes an error. Switching to relative rank.\n');
        rank_X = sum(abs(diag(D_X)) >= epsilon2*max(abs(diag(D_X))));
        if rank_X >= 3
            % if relative tolerance produces incorrect rank, then please
            % reselect tolerance epsilon2
            printf('The given epsilon2 produces unexpected rank of X* in both absolute and relative tolerance.\n Please reselect epsilon2.\n');
            status.exception.code = 1; % 1 for exception (rank(X_star) > 2) happens. Theoretically, rank(X_star) <= 2.
            status.exception.M0 = M0;
            status.exception.M1 = M1;
            status.exception.M2 = M2;
            status.exception.epsilon2 = epsilon2;
            return % return the control to the invoking function
        else
            useRelRankTolTag = 1;
        end
    end
    status.rank_X = rank_X;
    
    % compute the numerical rank of Z with tolerance epsilon2
    rank_Z = sum(abs(diag(D_Z)) >= epsilon2);
    status.rank_Z = rank_Z;
    
    % check condition 1
    if (y_hat1 > epsilon2) && (y_hat2 > epsilon2)
        if verbosity == 1
            fprintf('rank(Z*) = %d and rank(X*) = %d\n',rank_Z,rank_X);
        end
        
        % check condition 2
        if rank_Z == n-1
            % check conditions 3 and 4
            if rank_X == 2
                % conduct a rank-one decomposition of X*
                index = 1:(n+1);
                if ~useRelRankTolTag
                    index(diag(D_X) < epsilon2) = [];
                else
                    index(diag(D_X) < epsilon2*max(abs(diag(D_X)))) = [];
                end
                
                u1 = sqrt(D_X(index(1),index(1)))*V_X(:,index(1));
                u2 = sqrt(D_X(index(2),index(2)))*V_X(:,index(2));
                
                a = trace(M1*u1*u1');
                b = 2*trace(M1*u2*u1');
                c = trace(M1*u2*u2');
                
                if abs(a) >= epsilon2
                    t = (-b+sqrt(b^2-4*a*c))/2/a;
                    x_star1 = (t*u1+u2)/sqrt(t^2+1);
                    x_star2 = (-u1+t*u2)/sqrt(t^2+1);
                else
                    % in case u1 and u2 are the rank-one decomposition
                    % which are active at the first constraint
                    x_star1 = u1;
                    x_star2 = u2;
                end
                
                % store the value of a tentative rank-one decomposition
                status.SP_x1 = x_star1;
                status.SP_x2 = x_star2;
                
                % check the remaining conditions in 3 and 4
                if (trace(M1*x_star1*x_star1') < epsilon2) && ...
                        (trace(M1*x_star2*x_star2') < epsilon2) && ...
                        (abs(trace(M2*x_star1*x_star1')) > epsilon2) && ...
                        (abs(trace(M2*x_star2*x_star2')) > epsilon2) && ...
                        (abs(trace(M1*x_star1*x_star2')) > epsilon2)
                    
                    noOptimalityGapTag = 0;
                    if verbosity == 1
                        fprintf('An optimality gap exists.\n');
                    end
                    status.noOptimalityGap = 0;
                else
                    noOptimalityGapTag = 4;
                    if verbosity == 1
                        fprintf('There is no optimality gap.\n');
                    end
                    status.noOptimalityGap = 1;
                end
            else
                noOptimalityGapTag = 3;
                if verbosity == 1
                    fprintf('There is no optimality gap.\n');
                end
                status.noOptimalityGap = 1;
            end
        else
            noOptimalityGapTag = 2;
            if verbosity == 1
                fprintf('There is no optimality gap.\n');
            end
            status.noOptimalityGap = 1;
        end
    else
        noOptimalityGapTag = 1;
        if verbosity == 1
            fprintf('There is no optimality gap.\n');
        end
        status.noOptimalityGap = 1;
    end
    
    status.exception.code = 0; % 0 for no exception happening
    if noOptimalityGapTag ~= 0
        if rank_X == 1
            % conduct rank-one decomposition
            index = 1:(n+1);
            index(diag(D_X) < epsilon2) = [];
            
            x_star1 = sqrt(D_X(index(1),index(1)))*V_X(:,index(1));
            
            status.SP_x1 = x_star1;
        elseif rank_X == 2
            % conduct rank-one decomposition
            index = 1:(n+1);
            if ~useRelRankTolTag
                index(diag(D_X) < epsilon2) = [];
            else
                index(diag(D_X) < epsilon2*max(abs(diag(D_X)))) = [];
            end
            u1 = sqrt(D_X(index(1),index(1)))*V_X(:,index(1));
            u2 = sqrt(D_X(index(2),index(2)))*V_X(:,index(2));
            
            % check if any constraint is inactive
            if abs(y_hat1) < epsilon2 && abs(y_hat2) < epsilon2
                % both constraints inactive
                x_star1 = u1;
                x_star2 = u2;
            elseif abs(y_hat1) > epsilon2 && abs(y_hat2) < epsilon2
                % second constraint inactive
                % align with first constraint
                a = trace(M1*u1*u1');
                b = 2*trace(M1*u2*u1');
                c = trace(M1*u2*u2');
                
                t = (-b+sqrt(b^2-4*a*c))/2/a;
                x_star1 = (t*u1+u2)/sqrt(t^2+1);
                x_star2 = (-u1+t*u2)/sqrt(t^2+1);
            elseif abs(y_hat1) < epsilon2 && abs(y_hat2) > epsilon2
                % first constraint inactive
                % align with second constraint
                a = trace(M2*u1*u1');
                b = 2*trace(M2*u2*u1');
                c = trace(M2*u2*u2');
                
                t = (-b+sqrt(b^2-4*a*c))/2/a;
                x_star1 = (t*u1+u2)/sqrt(t^2+1);
                x_star2 = (-u1+t*u2)/sqrt(t^2+1);
            else
                % both constraints active
                % check if rankZ < n-1
                if rank_Z < n-1
                    % first check if u1 and u2 provide a pair of solution
                    if (abs(trace(M1*u1*u1')) < epsilon2) && ...
                            (abs(trace(M2*u1*u1')) < epsilon2) && ...
                            (abs(trace(M1*u2*u2')) < epsilon2) && ...
                            (abs(trace(M2*u2*u2')) < epsilon2)
                        
                        x_star1 = u1;
                        x_star2 = u2;
                    else
                        % try to align the rank-one decomposition with the first
                        % constraint
                        a = trace(M1*u1*u1');
                        b = 2*trace(M1*u2*u1');
                        c = trace(M1*u2*u2');
                        
                        t = (-b+sqrt(b^2-4*a*c))/2/a;
                        x_star1 = (t*u1+u2)/sqrt(t^2+1);
                        x_star2 = (-u1+t*u2)/sqrt(t^2+1);
                        
                        if  (abs(trace(M2*x_star1*x_star1')) <= epsilon2) && ...
                                (abs(trace(M2*x_star2*x_star2')) <= epsilon2)
                            % then such x_star1 and x_star2 can be a pair of primal
                            % optimal solution
                            % do nothing
                        elseif ((abs(trace(M2*x_star1*x_star1')) > epsilon2) || ...
                                (abs(trace(M2*x_star2*x_star2')) > epsilon2)) && ...
                                abs(trace(M1*x_star1*x_star2')) <= epsilon2
                            % M1*x_star1*x_star2' is zero
                            % re-align x_star1 and x_star2 with the second
                            % constraint
                            u1 = x_star1;
                            u2 = x_star2;
                            
                            a = trace(M2*u1*u1');
                            b = 2*trace(M2*u2*u1');
                            c = trace(M2*u2*u2');
                            
                            t = (-b+sqrt(b^2-4*a*c))/2/a;
                            x_star1 = (t*u1+u2)/sqrt(t^2+1);
                            x_star2 = (-u1+t*u2)/sqrt(t^2+1);
                        elseif ((abs(trace(M2*x_star1*x_star1')) > epsilon2) || ...
                                (abs(trace(M2*x_star2*x_star2')) > epsilon2)) && ...
                                abs(trace(M1*x_star1*x_star2')) > epsilon2
                            % in this case, a primal solution is generated
                            % from X_hat and Z_hat
                            
                            % sort out the non-zero eigenvalue associated
                            % eigenvectors of X_hat and Z_hat
                            % denoted by u_1,...,u_r,  r := rank(X_hat)
                            %            v_1,...,v_p,  p := rank(Z_hat)
                            % all column vectors
                            % note that r + p < n + 1
                            % [U,S,V] = svd([u_1,...,u_r,v_1,...,v_p]');
                            % then look for the columns of V such that the
                            % corresponding diag(S) is zero or the diagonal
                            % element is missing.
                            % These vectors lie in the intersection of
                            % N(X_hat) and N(Z_hat) (null space)
                            try
                                % collect the vectors that span the range space
                                % of X_hat
                                if ~useRelRankTolTag
                                    range_X = D_X(:,abs(diag(D_X)) >= epsilon2);
                                else
                                    range_X = D_X(:,abs(diag(D_X)) >= epsilon2*max(abs(diag(D_X))));
                                end
                                % collect the vectors that span the range space
                                % of Z_hat
                                range_Z = D_Z(:,abs(diag(D_Z)) >= epsilon2);
                                
                                % conduct a singular value decomposition to the
                                % vectors in range_X and range_Z
                                [U,S,V] = svd([range_X range_Z]');
                                
                                % collect the vectors that span the
                                % intersection of the null space of X_hat and
                                % null space of Z_hat.
                                null_X_Z = V(:,abs(diag(S)) < epsilon2);
                                
                                if size(S,1) < size(S,2)
                                    null_X_Z = [null_X_Z V(:,size(S,1)+1:size(S,2))];
                                end
                                
                                % take x_star3 that contains all directions in the null space
                                x_star3 = mean(null_X_Z,2);
                                
                                % then try to conduct a rank-one
                                % decomposition to x_star3*x_star3' + x_star1*x_star1'
                                % + x_star2*x_star2'
                                % to get a new x_star1 and x_star2
                                
                                % load the products
                                A133 = trace(M1*x_star3*x_star3');
                                A112 = trace(M1*x_star1*x_star2');
                                A113 = trace(M1*x_star1*x_star3');
                                A123 = trace(M1*x_star2*x_star3');
                                A211 = trace(M2*x_star1*x_star1');
                                A222 = trace(M2*x_star2*x_star2');
                                A233 = trace(M2*x_star3*x_star3');
                                A212 = trace(M2*x_star1*x_star2');
                                A213 = trace(M2*x_star1*x_star3');
                                A223 = trace(M2*x_star2*x_star3');
                                
                                % initialize parameters
                                alpha1 = [];
                                alpha2 = [];
                                alpha3 = [];
                                if abs(A112) < epsilon2
                                    % case 1
                                    alpha1 = 1;
                                    alpha3 = 0;
                                    alpha2 = (-2*A212+sqrt(4*A212^2-4*A222*A211))/2/A222;
                                else
                                    % case 2
                                    alpha3 = 1;
                                    t1 = -A123/A112;
                                    if abs(2*A113*A123-A112*A133) > epsilon2
                                        % case 2.1
                                        coefficient4 = 4*A112^2*A211;
                                        coefficient3 = 8*A112*A123*A211-8*A112*A113*A212+8*A112^2;
                                        coefficient2 = 4*A123^2*A211+4*A113^2*A222-4*A112*A133*A212-8*A123*A113*A212...
                                                       +16*A112*A123*A213-8*A113*A112*A223 + 4*A112^2*A233;
                                        coefficient1 = 4*A113*A133*A222-4*A123*A133*A212 + 8*A123^2*A213...
                                                       -8*A113*A123*A223-4*A133*A112*A223+8*A112*A123*A233;
                                        coefficient0 = A133^2*A222-4*A133*A123*A223+4*A123^2*A233;
                                        candidate_roots = roots([coefficient4 coefficient3 coefficient2 coefficient1 coefficient0]);
                                        for k = 1:length(candidate_roots)
                                            if isreal(candidate_roots(k)) && candidate_roots(k) >= t1
                                                % the solution required
                                                % needs to be real and also
                                                % greater or equal to t1
                                                alpha1 = candidate_roots(k);
                                                break;
                                            end
                                        end
                                        if isempty(alpha1)
                                            status.exception.code = 4; % 4 for alpha1 undetermined
                                            status.exception.M0 = M0;
                                            status.exception.M1 = M1;
                                            status.exception.M2 = M2;
                                            fprintf('Exception occurs. Please send the variable #status# to cheng@terpmail.umd.edu. Thank you!\n');
                                            return
                                        end
                                        alpha2 = -(2*alpha1*A113+A133)/2/(alpha1*A112+A123);  % p(alpha1)
                                    else
                                        % case 2.2
                                        a = A222^2;
                                        b = 2*t1*A212+2*A223;
                                        c = t1^2*A211 + 2*t1*A213 + A233;
                                        discriminant = b^2-4*a*c;
                                        if discriminant >= 0
                                            alpha2 = (-b + sqrt(discriminant))/2/a;
                                            alpha1 = t1;
                                        else
                                            t2 = -1/2;
                                            alpha2 = t2;
                                            a = A211;
                                            b = 2*t2*A212+2*A213;
                                            c = t2^2*A222+A233;
                                            discriminant = b^2-4*a*c;
                                            if (-b + sqrt(discriminant))/2/a > t1
                                                alpha1 = (-b + sqrt(discriminant))/2/a;
                                            else
                                                alpha1 = (-b - sqrt(discriminant))/2/a;
                                            end
                                        end
                                    end
                                end
                                x_star = alpha1*x_star1 + alpha2*x_star2 + alpha3*x_star3;
                                status.SP_x1 = x_star;
                            catch
                                % if there's an error, report it
                                status.exception.code = 3; % 3 for errors in the yet untested case
                                status.exception.M0 = M0;
                                status.exception.M1 = M1;
                                status.exception.M2 = M2;
                                fprintf('Exception occurs. Please send the variable #status# to cheng@terpmail.umd.edu. Thank you!\n');
                                return
                            end
                        end
                    end
                else
                    % first check if u1 and u2 provide a pair of solution
                    if (abs(trace(M1*u1*u1')) < epsilon2) && ...
                            (abs(trace(M2*u1*u1')) < epsilon2) && ...
                            (abs(trace(M1*u2*u2')) < epsilon2) && ...
                            (abs(trace(M2*u2*u2')) < epsilon2)
                        
                        x_star1 = u1;
                        x_star2 = u2;
                    else
                        % try to align the rank-one decomposition with the first
                        % constraint
                        a = trace(M1*u1*u1');
                        b = 2*trace(M1*u2*u1');
                        c = trace(M1*u2*u2');
                        
                        t = (-b+sqrt(b^2-4*a*c))/2/a;
                        x_star1 = (t*u1+u2)/sqrt(t^2+1);
                        x_star2 = (-u1+t*u2)/sqrt(t^2+1);
                        
                        if  (abs(trace(M2*x_star1*x_star1')) < epsilon2) && ...
                                (abs(trace(M2*x_star2*x_star2')) < epsilon2)
                            % then such x_star1 and x_star2 can be a pair of primal
                            % optimal solution
                            % do nothing
                        else
                            % M1*x_star1*x_star2' must be zero
                            % re-align x_star1 and x_star2 with the second
                            % constraint
                            u1 = x_star1;
                            u2 = x_star2;
                            
                            a = trace(M2*u1*u1');
                            b = 2*trace(M2*u2*u1');
                            c = trace(M2*u2*u2');
                            
                            t = (-b+sqrt(b^2-4*a*c))/2/a;
                            x_star1 = (t*u1+u2)/sqrt(t^2+1);
                            x_star2 = (-u1+t*u2)/sqrt(t^2+1);
                        end
                    end
                end
                
                if abs(abs(x_star1(1))-1) > abs(abs(x_star2(1))-1)
                    % choose the vector that has a relatively bigger first element
                    temp = x_star1;
                    x_star1 = x_star2;
                    x_star2 = temp;
                end
                
                status.SP_x1 = x_star1;
                status.SP_x2 = x_star2;
            end
        end
        
        % verification: solve the primal with initial value being the
        % rank-one decomposition result (using YALMIP)
        x_plus = sdpvar(n,1);
        options = sdpsettings('solver','fmincon','fmincon.TolFun',1e-9,'usex0',1,'verbose',verbosity);
        assign(x_plus,x_star1(2:end)/x_star1(1));
        optimize([[1;x_plus]'*M1*[1;x_plus]<=0;[1;x_plus]'*M2*[1;x_plus]<=0],[[1;x_plus]'*M0*[1;x_plus]],options);
        
        switch verbosity
            case 1
                fprintf('Primal minimum = %f, semidefinite relaxation optimum = %f.\n',value([1;x_plus]'*M0*[1;x_plus]),cvx_optval);
            case 0
                % do nothing
        end
        status.primal_optval = value([[1;x_plus]'*M0*[1;x_plus]]);  % optimal value of the primal problem
        status.primal_optsol = value(x_plus); % optimal solution of the primal problem
    end
end

