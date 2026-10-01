% Parameters
L = 15000;
dx = L/Nx;  % Set your desired dx value here
%n_lines = 100;  % Set the number of lines you want

start_at = 500;
%wavelength = 600;

end_at = start_at + wavelength*2;

n_lines = ceil(wavelength*2/dx)

% Open file for writing
  fid = fopen('/home/users/bpauken/suntans-gvc-2025/examples/quick-isw-test/rundata/dataxy.dat', 'w');

% Generate and write each line
for i = 1:n_lines
	  first_number = 500.0 + (i - 1) * dx;
second_number = 50.0;
if(n==n_lines)
  fprintf(fid, '%.1f %.1f', first_number, second_number);
 else
   fprintf(fid, '%.1f %.1f\n', first_number, second_number);
    end
end

% Close file
    fclose(fid);

disp('File generated successfully!');
