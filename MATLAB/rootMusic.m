function DOA_MuSiC =  rootMusic()
n = (0:99)';
frqs = [pi/4 pi/4+0.06];

r_0 = 2*exp(1j*frqs(1)*n)+1.5*exp(1j*frqs(2)*n)+ ...
    0.5*randn(100,1)+1j*0.5*randn(100,1);
        prediction_model_order = 12;
        [~,R_Music_0] = corrmtx(r_0,prediction_model_order,'modified');
        % W, a vector of frequencies in rad/sample.
        % P, corresponding signal power in the vector pow.
        subspace_dim = 2;
DOA_MuSiC = ang; str_rootMusic = ['DOA (Root-MuSiC) =  ' num2str(ang)];  disp(str_rootMusic);

end