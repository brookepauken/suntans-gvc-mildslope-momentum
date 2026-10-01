function omega = get_omega(dirname,n)

    addpath("~/suntans/suntans-master/mfiles/");

    suntansfile = [dirname,'/suntans.dat'];

    depthdata = load([dirname,'/depth.dat-voro'],'-ascii');

    dv = depthdata(:,3);

    Nkmax = getvalue(suntansfile,'Nkmax');
    Nc = length(dv);
    precision='float64';

    fid = fopen([dirname,'/vert_omega.dat'],'rb');
    %fseek(fid,8*(Nkmax+1)*Nc*(n-1),'bof');
    omega = reshape(getcdata(fid,Nc*(Nkmax+1),n,precision),Nc,Nkmax+1);
    %omega = reshape(fread(fid,(Nkmax+1)*Nc,'float64'),Nc,Nkmax+1);
    fclose(fid);
    % 
    % fid = fopen([dirname,'/fs.dat'],'rb');
    % fseek(fid,8*Nc*(n-1),'bof');
    % h = fread(fid,Nc,'float64');
    % fclose(fid);
    % 
    % omega = zeros(Nc,Nkmax+1);
    % omega(:,Nkmax+1) = 0;
    % %zeta(:, 1) = h*0;
    % for k=Nkmax:-1:1
    %     omega(:,k) = omega(:,k-1)-dzz(:,k-1);
    % end
end