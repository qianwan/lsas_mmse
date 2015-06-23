function UEs = brownian(L, K, BSs, outerRadius)
    UEs = zeros(L * K, 1);
    for q = 1 : L
        for k = 1 : K
            while true
                x = (rand - 0.5) * 2 * outerRadius;
                y = (rand - 0.5) * 2 * outerRadius;
                mx = abs(x);
                my = abs(y);
                valid = true;
                %if my > outerRadius * sin(pi / 3) || (sqrt(3) * mx + my > outerRadius * sqrt(3))
                if my > outerRadius * sin(pi / 3) || mx > outerRadius || (mx > outerRadius / 2 && my > (outerRadius - mx) * sin(pi / 3))
                     valid = false;
                end
                if valid == true
                    break;
                end
            end
            UEs((q - 1) * K + k) = x + y * 1j + BSs(q);
        end
    end
    return
