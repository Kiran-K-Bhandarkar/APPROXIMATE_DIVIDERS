% Main function
clc;
clear;

% Input parameters
n = 8;
a = 4567  ;     % Dividend
b = 162;            % Divisor
p = 6 ;              % p: Tuning parameter - p >= n/2

% Exact divider output
[qnt_ex, rem_ex] = exact_divider(dec2bin(a, 2*n), ...
    dec2bin(b, n), '0');

% In-exact divider output
[qnt_iex, rem_iex] = inexact_divider(dec2bin(a, 2*n), ...
    dec2bin(b, n), '0', p, n);

% Computing error distance between exact and inexact ouptuts (quotient)
Err_qnt = comparator(qnt_ex, qnt_iex);

% Computing error distance between exact and inexact ouptuts (reminder)
Err_rem = comparator(rem_ex, rem_iex);

% Display Output
fprintf('%i x %i Divider\n', 2*n, n);
fprintf('Dividend = %i, Divisor= %i, Tuning Factor = %i\n', a, b, p);
q_ex=bin2dec(qnt_ex);
r_ex=bin2dec(rem_ex);
fprintf('Exact Division:       Quotient= %i, Reminder= %i\n', ...
    q_ex, r_ex);
q_iex=bin2dec(qnt_iex);
r_iex=bin2dec(rem_iex);
fprintf('Approximate Division: Quotient= %i, Reminder= %i\n', ...
    q_iex, r_iex);
fprintf('Quotient Error distance: %i\n', (q_ex - q_iex));
fprintf('Reminder Error distance: %i\n', abs(r_ex - r_iex));

% ----------------- FUNCTIONS -------------------------%
% Main exact divider function
function [qnt, rem] = exact_divider(dividend, divisor, bin)
    l = length(divisor);
    for i = 1: l
        if i == 1
            rem_temp = dividend(1:(l+1));
        else
            rem_temp = strcat(r, dividend(l+i));
        end
        [q, r] = ex_divder_row(rem_temp, divisor, bin);
        qnt(i) = q;  %#ok<AGROW> 
        r = num2str(r);
        r = r(~isspace(r));   
%         disp(qnt(i));
%         disp(r);
    end
    qnt = num2str(qnt);
    qnt = qnt(~isspace(qnt));
    rem = r;   
end

% Main inexact divider function
function [qnt, rem] = inexact_divider(dividend, divisor, bin, p, n)
    l = length(divisor);
    for i = 1: l
        if i == 1
            rem_temp = dividend(1:(l+1));
        else
            rem_temp = strcat(r, dividend(l+i));
        end

        if (2*n-i)-p >= l
            [q, r] = ex_divder_row(rem_temp, divisor, bin);
        else
            if (2*n-i)-p > 0
                p_f = (2*n-i)-p;
            else
                p_f = 0;
            end
            [q, r] = iex_divder_row(rem_temp, divisor, bin, p_f);
        end
        qnt(i) = q;  %#ok<AGROW> 
        r = num2str(r);
        r = r(~isspace(r));   
    end
    qnt = num2str(qnt);
    qnt = qnt(~isspace(qnt));
    rem = r;   
end


% Function row for exact divider cell
function [qnt, rem] = ex_divder_row(dvd, div, bin)
    % Calculating the quotient bits
    l = length(div);
    b_temp = zeros(1, l);
    for i = 1: l
        if i == 1
            b_temp(i) = exrdc_b(bin2dec(dvd(l+1)), bin2dec(div(l)), ...
                bin2dec(bin));
        else
            b_temp(i) = exrdc_b(bin2dec(dvd(l-i+2)), bin2dec(div(l-i+1)), ...
                b_temp(i-1));
        end
    end
    qnt = or(bin2dec(dvd(1)),not(b_temp(l)));

    % Calculating the reminder bits
    rem = zeros(1, l);
    for i = 1: l
        if i == 1
            rem(l-i+1) = exrdc_r(bin2dec(dvd(l+1)), bin2dec(div(l)), ...
                bin2dec(bin), qnt);
        else 
            rem(l-i+1) = exrdc_r(bin2dec(dvd(l-i+2)), bin2dec(div(l-i+1)), ...
                b_temp(i-1), qnt);
        end
    end
end

% Function row for inexact divider cell
function [qnt, rem] = iex_divder_row(dvd, div, bin, p)
    % Calculating the quotient bits
    l = length(div);
    b_temp = zeros(1, l);
    for i = 1: l
        if i <= l-p
            if i == 1
                b_temp(i) = inexrdc_b(bin2dec(dvd(l+1)), bin2dec(div(l)), ...
                    bin2dec(bin));
            else
                b_temp(i) = inexrdc_b(bin2dec(dvd(l-i+2)), bin2dec(div(l-i+1)), ...
                    b_temp(i-1));
            end
        else
            if i == 1
                b_temp(i) = exrdc_b(bin2dec(dvd(l+1)), bin2dec(div(l)), ...
                    bin2dec(bin));
            else
                b_temp(i) = exrdc_b(bin2dec(dvd(l-i+2)), bin2dec(div(l-i+1)), ...
                    b_temp(i-1));
            end
        end

    end
    qnt = or(bin2dec(dvd(1)),not(b_temp(l)));

    % Calculating the reminder bits
    rem = zeros(1, l);
    for i = 1: l
        if i <= l-p
            if i == 1
                rem(l-i+1) = inexrdc_r(bin2dec(dvd(l+1)), bin2dec(div(l)), ...
                    bin2dec(bin), qnt);
            else 
                rem(l-i+1) = inexrdc_r(bin2dec(dvd(l-i+2)), bin2dec(div(l-i+1)), ...
                    b_temp(i-1), qnt);
            end
        else
            if i == 1
                rem(l-i+1) = exrdc_r(bin2dec(dvd(l+1)), bin2dec(div(l)), ...
                    bin2dec(bin), qnt);
            else 
                rem(l-i+1) = exrdc_r(bin2dec(dvd(l-i+2)), bin2dec(div(l-i+1)), ...
                    b_temp(i-1), qnt);
            end
        end
    end
end


% Function for exact divder cell (borrow) - basic block of the divder array
function bout = exrdc_b(x, y, bin)
    term1 = and(not(xor(x,y)),bin);
    term2 = and(not(x),y);
    bout = or(term1, term2); % Borrow term in subtractor cell
end

% Function for exact divder cell (reminder) - basic block of divder array
function r = exrdc_r(x, y, bin, qs)
    d = xor(xor(x, y), bin); % Difference term in subtractor cell
    if qs == 1
        r = d;
    else
        r = x;
    end
end

% Function for inexact divder cell (borrow) - basic block of divder array
function bout = inexrdc_b(x, y, bin) %#ok<INUSD> 
    bout = not(x);
end

% Function for inexact divder cell (reminder) - basic block of divder array
function r = inexrdc_r(x, y, bin, qs) %#ok<INUSL> 
    r = xor(qs, x);
end

% Function to compared the two inputs and output the distance
function Err = comparator(A, B)
    A_dec = bin2dec(A);
    B_dec = bin2dec(B);
    Err = abs(A_dec-B_dec);
end

