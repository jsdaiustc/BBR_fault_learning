function [mu,Pm,f_sample_all]=VBI_offgrid_CGDP(Y,C_all,f_sample,Fs)
Y = reshape( Y,[length(Y),1] );
[M, Snap] = size(Y);
norm_y = norm(Y,'fro') / sqrt(M*Snap);
Y = Y / norm_y;
f_sample_all = f_sample;
%% Initialization
A = exp( 1i*2*pi*C_all/Fs*f_sample_all ) / sqrt(M);
B = diag( 1i*2*pi*C_all/Fs ) * A;
reslu = f_sample(2) - f_sample(1);
N = size(f_sample,2);
converged = false;
iter = 0;
mu = 0;
delta = ones(N,1);
maxiter = 450;
alpha0 = 1;
etc = 10;
xi = ones(N,1);
h = 0.001;
delta_inv = 1./delta;
while ~converged
   %% update mu and sigma
    AA = A' * A;
    Sigma = inv( alpha0* AA  + diag(delta) );
    mu = Sigma * ( alpha0 * (A' * Y) );    
   %% update delta
    delta = 1 ./ delta_inv;
    sigma_diag = diag(Sigma);
    z = (N+1) / 2;
    part1 = abs(mu(z+1:end)).^2;
    part2 = abs(mu(z-1:-1:1)).^2;
    part3 = sigma_diag(z+1:end);
    part4 = sigma_diag(z-1:-1:1);
    temp1 = part1 + part2 + part3 + part4;
    temp2 = abs(mu(z))^2 + sigma_diag(z);
    delta_inv(z+1:end) = ( -3+2*sqrt(9/4+xi(z+1:end).^2.*temp1) )./( xi(z+1:end).^2 );
    delta_inv(z-1:-1:1) = delta_inv(z+1:end);
    delta_inv(z) = ( -1+sqrt(1+4*xi(z)^2*temp2) ) / (xi(z)^2);
   %% update xi
    xi((N+1)/2:end) = ( -h+sqrt( h^2+2*delta_inv((N+1)/2:end)*(h+2) ) )./(delta_inv((N+1)/2:end));
    xi(1:(N+1)/2-1) = fliplr( xi((N+1)/2+1:end)' );   
   %% update beta
    resid = Y - A*mu;
    beta_old = alpha0;   
    alpha0 = (M*Snap) / (((norm(Y-A*mu,'fro'))^2) + Snap*real( sum(diag( Sigma*AA)) ) );
    rho_beta = 0.95;
    alpha0 = (1-rho_beta)*alpha0 + rho_beta * beta_old ;    
    Pm = sum( mu.*conj(mu), 2 );
%     figure(5);stem(f_sample_all,10*log10(abs(Pm))+80 )
   %% grid refine
    idx = [1:N];
    BHB = B(:,idx)' * B(:,idx);
    P = real(conj(BHB) .* (mu(idx,:) * mu(idx,:)' + Snap * Sigma(idx,idx)));
    v = sum( real(conj(mu(idx,:)) .* (B(:,idx)' * (resid))),2);
    v = v - Snap * real(diag(B(:,idx)' * A * Sigma(:,idx)));
    temp_grid = sum(v) / sum(sum(P));
    increa_grid = sign(temp_grid) / 1000 * 0.95^(iter);
    f_sample_all(idx) = f_sample_all(idx) + increa_grid;
    A(:,idx) = exp( 1i*2*pi*C_all / Fs * f_sample_all(idx) ) / sqrt(M);
    B(:,idx) = diag( 1i*2*pi*C_all/Fs ) * A(:,idx);
    choose_set = 1:N;
    Pm = sum( mu(choose_set,:) .* conj(mu(choose_set,:)), 2);
    Pm = Pm(1:(end-1)/2)  +  Pm(end:-1:(end-1)/2+2);
    [~,sort_ind] = sort(Pm, 'descend');
    etc = min(length(Pm),etc);
    idx = sort_ind(1:etc);
    idx = [idx;(N+1)-idx];
    BHB = B(:,idx)' * B(:,idx);
    P = real(conj(BHB) .* (mu(idx,:) * mu(idx,:)' + Snap * Sigma(idx,idx)));
    v = sum( real(conj(mu(idx,:)) .* (B(:,idx)' * (Y-A*mu))),2);
    v = v - Snap * real(diag(B(:,idx)' * A * Sigma(:,idx)));
    Permatrix = [-eye(etc);eye(etc) ];
    P = Permatrix'*P*Permatrix;
    v = Permatrix'*v;
    temp_grid = v./diag(P);
    temp_grid = temp_grid';
    theld = reslu * 0.95^(iter);
    ind_small = find( abs(temp_grid)<theld );
    temp_grid(ind_small) = sign( temp_grid(ind_small) ) * theld;
    ind_unchang = find( abs(temp_grid)>reslu );
    temp_grid(ind_unchang) = sign( temp_grid(ind_unchang) ) * reslu/20;
    f_sample_all(idx) = f_sample_all(idx) + [-temp_grid(1:etc), temp_grid(1:etc)];
    A(:,idx) = exp( 1i*2*pi*C_all/Fs*f_sample_all(idx) ) /sqrt(M);
    B(:,idx) = diag(1i*2*pi*C_all/Fs)*A(:,idx);
    %% stopping criteria
    if iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
end
Pm = sum( mu.*conj(mu), 2);
[f_sample_all, indsort] = sort(f_sample_all);
Pm = Pm(indsort);







