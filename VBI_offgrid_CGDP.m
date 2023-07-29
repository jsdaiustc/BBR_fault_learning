function [mu,Pm,f_sample]=VBI_offgrid_CGDP(y,C_all,f_sample,Fs)
y = reshape( y,[length(y),1] );
M = length(y);
N = length(f_sample);
norm_y = norm(y,'fro') / sqrt(M);
y = y / norm_y;
%% Initialization
omega = 1i*2*pi * C_all / Fs;
A = exp( omega * f_sample ) / sqrt(M);
C = 1i*2*pi * C_all / Fs * ones( 1, N );
B = C .* A;
reslu = f_sample(2) - f_sample(1);
converged = false;
iter = 0;
mu = 0;
delta = ones(N,1);
maxiter = 150;
alpha = 1;
etc = 10;
xi = ones(N,1);
h = 2e-4;
a = 1e-5;
b = a;
delta_inv = 1./delta;
while ~converged
   %% update mu and sigma
    AHA = A' * A;
    Sigma = ( alpha* AHA  + diag(delta) ) \ eye(N);
    mu = Sigma * ( alpha * (A' * y) );    
   %% update delta
    delta = 1 ./ delta_inv;
    sigma_diag = diag(Sigma);
    z = (N+1) / 2;
    part1 = abs( mu(z+1: end) ).^2;
    part2 = abs( mu(z-1: -1: 1) ).^2;
    part3 = sigma_diag( z+1: end );
    part4 = sigma_diag( z-1: -1: 1 );
    temp1 = part1 + part2 + part3 + part4;
    temp2 = abs(mu(z))^2 + sigma_diag(z);
    delta_inv(z+1:end) = ( -3+2*sqrt( 9/4 + xi(z+1:end).^2.*temp1 ) ) ./ ( xi(z+1:end).^2 );
    delta_inv(z-1:-1:1) = delta_inv( z+1: end );
    delta_inv(z) = ( -1 + sqrt( 1 + 4 * xi(z)^2 * temp2 ) ) / (xi(z)^2);
   %% update xi
    xi((N+1)/2:end) = ( -h + sqrt( h^2 + 2*delta_inv((N+1)/2:end)*(h+2) ) ) ./ (delta_inv((N+1)/2:end));
    xi(1:(N+1)/2-1) = fliplr( xi( (N+1) / 2 + 1: end )' );   
   %% update alpha
    resid = y - A * mu;
    alpha_old = alpha;   
    alpha = ( a + M ) / ( b + norm(resid,'fro')^2 + real( sum( sum(conj(AHA).*Sigma) ) ) );   
    rho_alpha = 0.95;
    alpha = ( 1 - rho_alpha ) * alpha + rho_alpha * alpha_old ;    
   %% grid refine
    BHB = B' * B;
    varpi = mu * mu' + Sigma;
    P = conj(BHB) .* varpi;
    P = real(P);
    v = real(  conj( mu ) .* ( B' * resid )  );
    v = v - real( diag(B' * A * Sigma) );    
    temp_grid = sum(v) / sum(sum(P));
    increa_grid = sign(temp_grid) / 1000 * 0.95^(iter);
    f_sample = f_sample + increa_grid;
    A = A .* exp( omega * increa_grid );
    B = C .* A;
    Pm = sum( mu .* conj(mu), 2);
    Pm = Pm( 1: (end-1) / 2 ) + Pm( end: -1: (end-1) / 2 + 2 );
    [~, sort_ind] = sort(Pm, 'descend');
    idx = sort_ind(1:etc);
    idx = [ idx; (N+1) - idx ];
    BHB = B(:,idx)' * B(:,idx);
    P = conj(BHB) .* varpi(idx,idx);  
    P = real(P);    
    v = real( conj( mu(idx,:) ) .* ( B(:,idx)' * resid ) );
    v = v - real( diag( B(:,idx)' * A * Sigma(:,idx) ) );
    Permatrix = [-eye(etc); eye(etc) ];
    P = Permatrix' * P * Permatrix;
    v = Permatrix' * v;
    temp_grid = v ./ diag(P);
    temp_grid = temp_grid';
    theld = reslu * 0.95^(iter);
    ind_small = find( abs( temp_grid ) < theld );
    temp_grid(ind_small) = sign( temp_grid(ind_small) ) * theld;
    ind_unchang = find( abs(temp_grid)>reslu );
    temp_grid(ind_unchang) = sign( temp_grid(ind_unchang) ) * reslu / 20;
    f_sample(idx) = f_sample(idx) + [-temp_grid(1:etc), temp_grid(1:etc)];
    A(:,idx) = A(:,idx) .* exp( omega *[-temp_grid(1:etc), temp_grid(1:etc)] );
    B(:,idx) = C(:,1:length(idx)) .* A(:,idx);
    %% stopping criteria
    if iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
end
Pm = sum( mu.*conj(mu), 2);
[f_sample, indsort] = sort(f_sample);
Pm = Pm(indsort);