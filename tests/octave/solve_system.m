% solve_system.m
% Usage: octave --quiet tests/octave/solve_system.m matrix.txt vector.txt out.txt
args = argv();
if numel(args) < 3
  error('Usage: solve_system.m matrix_file vector_file out_file');
end
matrix_file = args{1};
vector_file = args{2};
out_file = args{3};
A = dlmread(matrix_file);
b = dlmread(vector_file);
x = A \ b;
dlmwrite(out_file, x, 'delimiter', ' ');
exit;
