%% Problem 2 part (i)
space = 0.5;
alpha = 2;
alpha2 = alpha^2;
K = 225;
n = (sqrt(K) - 1) / 2;

% enumerating the bases
i = 1;
k_sqr = [0, 0];
for y = n*space:-space:-n*space
    for x = -n*space:space:n*space
        if x == y
            break;
        end
        k_sqr(i, :) = [x, y];
        i = i + 1;
    end
end

% calculating the matrix elements
num_pts = length(k_sqr);
S = zeros(num_pts, num_pts);
Hx = zeros(num_pts, num_pts);
Hy = zeros(num_pts, num_pts);
A = zeros(num_pts, num_pts);
for ia = 1:num_pts
    xA = k_sqr(ia, 1);
    yA = k_sqr(ia, 2);
    for ib = 1:num_pts
        xB = k_sqr(ib, 1);
        yB = k_sqr(ib, 2);

        s_xaxb = s(xA,xB, alpha);
        s_yayb = s(yA,yB, alpha);
        s_xayb = s(xA,yB, alpha);
        s_yaxb = s(yA, xB, alpha);
        S(ia, ib) = 2*s_xaxb*s_yayb - 2*s_xayb*s_yaxb;
        
        f_xaxb = f(xA,xB, alpha);
        f_yayb = f(yA,yB, alpha);
        f_xayb = f(xA,yB, alpha);
        f_yaxb = f(yA, xB, alpha);
        Hx(ia, ib) = f_xaxb*s_yayb - f_yaxb*s_xayb ...
             - f_xayb*s_yaxb + f_yayb*s_xaxb;
        Hy(ia, ib) = s_xaxb*f_yayb - s_yaxb*f_xayb ...
                     - s_xayb*f_yaxb + s_yayb*f_xaxb;
                 
        A(ia, ib) = 2*g(xA,xB,alpha)*g(yA,yB,alpha) ...
                    - 2*g(xA,yB,alpha)*g(yA,xB,alpha);
        
    end
end

%% Problem 2 part (ii)
echodemo('prob2.m', 1);

fprintf('Exact Ground State Energies:\n');
aList = [0, 0.001, 0.01, 0.1, 1, 10];
exact_E0 = zeros(1, length(aList));
for i = 1:length(aList)
    a = aList(i);
    [wf, D] = eig(S \ (Hx + Hy + a*A));
    energ = diag(D);
    [energ_sort, indices] = sort(energ);
    exact_E0(i) = energ_sort(1);
    
    fprintf('%f & %f \\\\ \n', aList(i), exact_E0(i));
end

%% Problem 2 part (iv)
echodemo('prob2.m', 1);

aList_p = 0:0.001:0.01;
Pert_E0 = zeros(1, length(aList));
for i = 1:length(aList_p)
    Pert_E0(i) = 2 + aList_p(i)*90/32;
end

aList_e = 0:0.001:0.01;
exact_E0 = zeros(1, length(aList_e));
for i = 1:length(aList_e)
    a = aList_e(i);
    [wf, D] = eig(S \ (Hx + Hy + a*A));
    energ = diag(D);
    [energ_sort, indices] = sort(energ);
    exact_E0(i) = energ_sort(1);
end

figure;
hold on;
plot(aList_e, exact_E0, 'o');
plot(aList_p, Pert_E0);
title('Anharmonic Oscillator with Antisymmetry');
xlabel('a');
xlabel('Energy');
legend('exact', 'perturbation');

%% Problem 2 part (v)
echodemo('prob2.m', 1);

% Normalizing basis coefficients
a = 10;
[wf, D] = eig(S \ (Hx + Hy + a*A));
energ = diag(D);
[energ_sort, indices] = sort(energ);

c = wf(:, indices(1));
norm = 1/sqrt(c'*S*c);
c = c * norm;

xvals = -2:0.1:2;
yvals = xvals;
psi = zeros(length(xvals), length(yvals));
for ix = 1:length(xvals)
    x = xvals(ix);
    for iy = 1:length(yvals)
        y = yvals(iy);
        for ia = 1:num_pts
            xA = k_sqr(ia, 1);
            yA = k_sqr(ia, 2);
            psi(ix,iy) = psi(ix,iy) + c(ia) * (exp(-alpha*(x-xA)^2 - alpha*(y-yA)^2) ...
                         - exp(-alpha*(x-yA)^2 - alpha*(y-xA)^2));
        end
    end
end

contourm(xvals, yvals, psi, 20);
xlabel('x');
ylabel('y');
title('Exact Anharmonic Oscillator Ground State Wavefn');