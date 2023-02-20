clear;
close;

n = 8;
Err_qnt_tot = [];
MRED_qnt_tot = [];

for t = 2: (2^n-1)
    Err_qnt = [];
    MRED_qnt = [];

    for w = 1: ((t*(2^n-1))+(2^n-1))
        a = w;           % Dividend
        b = t;           % Divisor
        r = 8 ;          % p: parameter (# pieces), can be 2,4,8
        R = 5 ;          % Rounding Width: 4,5...
    
        % Exact divider output
        qnt_ex = exact_div(a,b);
        % In exact divider output
        qnt_iex = inexact_div(a, b, r, R, n);
        
        % Computing error distance between exact and inexact ouptuts (quotient)
        Err_qnt(w) = comparator(floor(qnt_ex), floor(qnt_iex));
        MRED_qnt(w) = Err_qnt(w)/qnt_ex; 

    end
    Err_qnt_tot = [Err_qnt_tot, Err_qnt];
    MRED_qnt_tot = [MRED_qnt_tot, MRED_qnt];
 
end

% Metrics Calculation
Max_Err_Qnt = max(Err_qnt_tot);
Norm_Err_qnt = mean(Err_qnt_tot/Max_Err_Qnt);
MRED_qnt = mean(MRED_qnt_tot);
ER_qnt = 100 * (nnz(Err_qnt_tot)/numel(Err_qnt_tot));


%------------- FUNCTIONS ----------------------
function qnt = exact_div(a, b)
    [Fa, ka] = log2(a);
    [Fb, kb] = log2(b);
    Xb = Fb^-1;
    qnt = (Fa*Xb)*(2^(ka-kb));
end

function qnt = inexact_div(a, b, r, R, n)
    [Fa, ka] = log2(a); 
    [~, kb] = log2(b); 
    b_bin = dec2bin(b,n);
    fb = mux(b_bin, kb, n);
    frb = rnd(fb, R);
    Xb = inverse(frb, r);
    qnt = (Fa*Xb)*(2^(ka-kb+1));
end

function out = inverse(b, r)
    decimal_val = b;
    if r == 2
        if decimal_val >= 0 && decimal_val < 0.4961
            alpha = -0.6684;
            beta = 1;
        else
            alpha = -0.3349;
            beta = 0.8345;
        end
    elseif r == 4
        if decimal_val >= 0 && decimal_val < 0.2461
            alpha = -0.8025;
            beta = 1;
        elseif decimal_val >= 0.2461 && decimal_val < 0.4961
            alpha = -0.5364;
            beta = 0.9345;
        elseif decimal_val >= 0.4961 && decimal_val < 0.7461
            alpha = -0.3828;
            beta = 0.8583;
        else
            alpha = -0.2869;
            beta = 0.7868;
        end
    else
        if decimal_val >= 0 && decimal_val < 0.1211
            alpha = -0.892;
            beta = 1;
        elseif decimal_val >= 0.1211 && decimal_val <  0.2461
            alpha = -0.7158;
            beta = 0.9787;
        elseif decimal_val >= 0.2461 && decimal_val < 0.3711
            alpha = -0.5853;
            beta = 0.9465;
        elseif decimal_val >= 0.3711 && decimal_val < 0.4961
            alpha = -0.4875;
            beta = 0.9103;
        elseif decimal_val >= 0.4961 && decimal_val < 0.6211
            alpha = -0.4123;
            beta = 0.8730;
        elseif decimal_val >= 0.6211 && decimal_val < 0.7461
            alpha = -0.3533;
            beta = 0.8363;
        elseif decimal_val >= 0.7461 && decimal_val < 0.8711
            alpha = -0.3061;
            beta = 0.8011;
        else
            alpha = -0.2677;
            beta = 0.7677;
        end
    end
    out = (decimal_val * alpha) + beta;
end

% Multiplexer Block
function fb = mux(in_bin, idx, n)
    id = n + 1 - idx;
    temp = in_bin(id+1:end);
    zero_pad = num2str(zeros(1, (n - length(temp))));
    fb = strcat(temp, zero_pad(~isspace(zero_pad)));
end

% Rounding Unit
function fr = rnd(in_bin, r)
    temp = num2str(zeros(1, r-1));
    temp = strcat(temp(~isspace(temp)), '1');
    if in_bin(r+1) == '1'
        fr = dec2bin(bin2dec(in_bin(1:r))+bin2dec(temp), r);
        if length(fr) > r
            fr = fr(2:end);
        end
    else
        fr = in_bin(1:r);
    end
    zero_pad = num2str(zeros(1,length(in_bin) - r));
    fr_bin = strcat(fr, zero_pad(~isspace(zero_pad)));
    fr = (str2double(fr_bin(1))/2) + (str2double(fr_bin(2))/4) + ...
        (str2double(fr_bin(3))/8) + (str2double(fr_bin(4))/16) + ...
        (str2double(fr_bin(5))/32) + (str2double(fr_bin(6))/64) + ...
        (str2double(fr_bin(7))/128) + (str2double(fr_bin(8))/256);
end

% Function to compared the two inputs and output the distance
function Err = comparator(A, B)
    Err = abs(A-B);
end