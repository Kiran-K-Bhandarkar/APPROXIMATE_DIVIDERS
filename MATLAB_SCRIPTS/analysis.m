ssim_result = readmatrix("results_generated.txt", 'Delimiter','\t');
hsim_result = readmatrix("divider_output.txt", 'Delimiter','\t');

for i= 1: 9983
    % quotient check
    if ssim_result(i, 1) ~= hsim_result(i, 1)
        fprintf("Quotient mismatch at line %d, SW_Sim: %d, HW_Sim: %d\n", ...
            i, ssim_result(i, 1), hsim_result(i, 1))
    end
    % reminder check
    if ssim_result(i, 2) ~= hsim_result(i, 2)
        fprintf("Reminder mismatch at line %d, SW_Sim: %d, HW_Sim: %d\n", ...
            i, ssim_result(i, 2), hsim_result(i, 2))
    end
end


