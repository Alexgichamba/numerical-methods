%Jacobi iterative method
format long
%class example 1 (remove comment symbols to solve)
%aMatrix = [0 -2 10 ; 10 -1 0; -1 10 -2];
%bMatrix = [6 ; 9 ; 7];

%example 2 (remove comment symbols to solve)
%aMatrix = [-1 11 -1 3;2 -1 10 -1;0 -3 -1 8; 10 -1 2 0 ];
%bMatrix = [25; -11; 15; 6];

%example 3 (remove comment symbols to solve)
%aMatrix = [2 -1 1; 1 1 1; -1 -1 2];
%bMatrix = [-1 ; 2 ; -5];

%example 4 (remove comment symbols to solve)
%aMatrix = [3 6 2; 3 -1 1; 3 3 7];
%bMatrix = [0; 1; 4];

%checks for dimension errors
aSize = size(aMatrix);
bSize = size(bMatrix);
if (aSize(1) ~= aSize(2)) && (aSize(2)~= bSize(1))
    disp('Check matrix dimensions')
end

n = bSize(1);

%function that pivots the matrix
%copy matrices will be used during the pivoting
cpyLhs = aMatrix;
cpyRhs = bMatrix;
for a = 1:n
    [maxA,index] = max(cpyLhs(:,a));
    %swap rows index and 1 for both rhs and lhs
    cpyLhs([a index],:) = cpyLhs([index a],:);
    cpyRhs([a index],:) = cpyRhs([index a],:);
end
aMatrix = cpyLhs;
bMatrix = cpyRhs;

%Now the pivoted matrix can be solved by Jacobi's iterative method
%this loop uses the general equation for x(k+1)
max_iterations = 100;
x = zeros(n,1);
x_plus1 = zeros(n,1);
tolerance = 1*10^-8; %tolerance can be adjusted
for iteration = 1:max_iterations
    %set default value of convergence
    convergence = true;
    for i = 1:n
        sum = 0;
        for j = 1:n
            if(j~=i) %the general eqn doesnt apply for j=i
                sum = sum + aMatrix(i,j)*x(j);
            end
        end
       x_plus1(i) = (-1/aMatrix(i,i)) * (sum - bMatrix(i));
       if(abs(x_plus1-x)>tolerance)
           convergence = false;
       end
    end
    if convergence
        break
    end
    x = x_plus1;
end
disp('iteration count')
disp(iteration)
disp('solution')
disp(x_plus1)