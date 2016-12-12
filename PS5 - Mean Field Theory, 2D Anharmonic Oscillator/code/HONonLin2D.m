function [energ, wf] = HONonLin2D(K, alpha, space, anHarm)

    alpha2 = alpha^2;
    n = (sqrt(K) - 1) / 2;

    % enumerating the bases
    k = -n*space:space:n*space;
    k_sqr = zeros(length(k)^2, 2);
    i = 1;
    for x = -n*space:space:n*space
        for y = -n*space:space:n*space
            k_sqr(i, :) = [x, y];
            i = i + 1;
        end
    end

    % calculating the matrix elements
    num_pts = length(k_sqr);
    S = zeros(num_pts, num_pts);
    H = zeros(num_pts, num_pts);
    A = zeros(num_pts, num_pts);
    for ia = 1:num_pts
        xA = k_sqr(ia, 1);
        yA = k_sqr(ia, 2);
        for ib = 1:num_pts
            xB = k_sqr(ib, 1);
            yB = k_sqr(ib, 2);

            x_diff2 = (xA - xB)^2;
            y_diff2 = (yA - yB)^2;
            xP2 = ((xA + xB) / 2)^2;
            yP2 = ((yA + yB) / 2)^2;

            S(ia, ib) = (pi/ 2*alpha) * exp(-alpha*x_diff2/2 - alpha*y_diff2/2);
            H(ia, ib) = (S(ia,ib)/2) * ...
                (2*alpha - alpha2*(x_diff2 + y_diff2) + 1/(2*alpha) + xP2 + yP2);
            
            % anharmonic term
            A(ia, ib) = (3/(16*alpha2) + 3*xP2/(2*alpha) +xP2^2) * ...
                        (3/(16*alpha2) + 3*yP2/(2*alpha) +yP2^2) * ...
                        S(ia, ib);
        end
    end

    % solving for eigenvalues/eigenfns
    [wf, D] = eig(S \ (H + (anHarm * A)));
    energ = diag(D);
    for i = 1:size(wf, 2)
        c = wf(:, i);
        norm = 1/sqrt(c'*S*c);
        wf(:, i) = c * norm;
    end

end

