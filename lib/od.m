function y = od(H)
a = det(H' * H);
prod = 1;
for c = 1 : size(H, 2)
    prod = prod * (H(:, c)' * H(:, c));
end
y = 1 - a / prod;