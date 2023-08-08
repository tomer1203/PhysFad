function sum_rate = getSumRate(H,P,sigmaN)
    H_sz = size(squeeze(H(1,:,:)));
    freq_rate = zeros(1,length(H));
    for f=1:length(H)
        Hf = squeeze(H(f,:,:));
        [~,S,~] = svd(Hf);
        freq_rate(f) = log2(det(eye(H_sz(1)) + S*S'.*P(f)/sigmaN));
    end
    sum_rate = abs(sum(freq_rate));
end