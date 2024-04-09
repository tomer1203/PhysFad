function [avg_diversity,diversitySampleEigenVec,diversitySampleVec] = getDataDiversity(sampled_Hs)
     H_size = size(sampled_Hs);
     NumofEigen = min(H_size(3),H_size(4));
     eigenMatrix = zeros(H_size(1),H_size(2),NumofEigen);
     for sample = 1:H_size(1) 
         for freq = 1:H_size(2)
            [~,S,~] = svd(squeeze(sampled_Hs(sample,freq,:,:)));
            eigenMatrix(sample,freq,:) = diag(S);
         end
     end
     diversitySampleEigenVec = std(eigenMatrix,0,2); % calculate standard divation on the frequency
     diversityFreqEigenVec = std(eigenMatrix,0,1); % calculate standard divation on the samples
     diversitySampleVec = mean(diversitySampleEigenVec,3);
     avg_diversity = mean(diversityFreqEigenVec,"all");
end

