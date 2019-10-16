function [KSp, opts] = simulateUnmixing (optsin)
%This function generates uncorrelated latent response distributions and
%simulates the process of imaging, followed by unmixing and generation of spatially
%shuffled nulls from the imaging data.
%This function was used to validate null distributions for our two-color imaging analyses

%output:    KSp: the p-value of a two-sample KS test between the null and the measured Z-scored correlation distribution
%           opts: the user-selected options

%Kaspar Podgorski 2019

opts.unmixHow = {'''LS''', '''RL'''}; tooltips.unmixHow = 'Least Squares or Richardson Lucy unmixing?';
opts.bA = 5; tooltips.bA = 'brightness of indicator A, ~photons/trial';
opts.bB = 10; tooltips.bB = 'brightness of indicator B, ~photons/trial';
opts.dist = {'''Uniform''', '''Exponential''', '''Gaussian''', '''Spike/Slab'''}; tooltips.dist = 'Simulated activity distribution';
opts.mix = [0.08 0.11]; tooltips.mix =  'mixing coefficients [A->B, B->A]';
opts.addGlobal = true; tooltips.addGlobal = 'Add a global signal to increase all correlations';
opts.randSeed = 42; tooltips.randSeed = 'For repeatable results, you can specify the pseudorandom seed';

if nargin %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
else
    opts = optionsGUI(opts, tooltips);
end

rng(opts.randSeed);

N = 1e4;
Ntrials = 200;

%simulate uncorrelated ("null-distributed") latents; 'A' and 'B' channels
switch opts.dist
    case 'Uniform'
        A = 2*rand(N, Ntrials); %synapses, trials
        B = 2*rand(N, Ntrials); %synapses, trials
    case 'Exponential'
        A = 2*rand(N, Ntrials); %synapses, trials
        B = 2*rand(N, Ntrials); %synapses, trials
        A = exp(A)-1;
        B = exp(B)-1;
    case 'Gaussian'
        A = max(0, 3+randn(N, Ntrials))/3; %synapses, trials
        B = max(0, 3+randn(N, Ntrials))/3; %synapses, trials
    case 'Spike/Slab'
        A = 3*rand(N, Ntrials)-1; %synapses, trials
        B = 3*rand(N, Ntrials)-1; %synapses, trials
        A(A<0) = 0;
        B(B<0) = 0;
end

%add global signals (shift correlations)
if opts.addGlobal
    G = rand(1, Ntrials)/10;
    A = A+G;
    B = B+G;
end

A = opts.bA.*A;
B =opts.bB.*B;
GT = [A(:) B(:)];

mix = opts.mix;
MMgt = [1-mix(1),mix(1);mix(2),1-mix(2)];%ground truth mixing matrix
C = poissrnd(GT*MMgt);

%mixing matrix used for recovery
unmix = mix.*(rand(size(mix))+0.5); %simulate unmixing matrix errors of up to 50%
MM = [1-unmix(1),unmix(1);unmix(2),1-unmix(2)];

y = C;
switch opts.unmixHow
    case 'RL' %unmix RL
        numIters = 40;
        darkRate = 1e-2;
        U = ones(size(y)); %initialize unmixed
        
        YP = repmat(sum(MM,2)', size(y,1),1);
        for iter = 1:numIters
            CC = (U+darkRate)*MM; %forward model of guess
            YC = (max(0,(y./CC)))*MM'; %unmix of normalized data
            U = U.*YC./YP;
            U(U<0) = 0; U(isnan(U)) = 0;
        end
        
        T = U';
    otherwise %unmix least squares
        T = MM'\y';
end

%Estimates of A and B
Au = reshape(T(1,:), size(A));
Bu = reshape(T(2,:), size(A));

%Compute sample correlations
Au = Au-mean(Au,2);
Bu = Bu-mean(Bu,2);
CC = sum(Au.*Bu,2)./(sqrt(sum(Au.^2,2)).*sqrt(sum(Bu.^2,2))); %sample correlation

%Compute shuffled null
tic;
disp('Computing correlations for 1000 shuffles...')
for shuff = 1000:-1:1
    if ~mod(shuff, 100)
        disp(shuff)
    end
    ixs = randperm(size(Au,1));
    CCss(:,shuff) = sum(Au(ixs,:).*Bu,2)./(sqrt(sum(Au(ixs,:).^2,2)).*sqrt(sum(Bu.^2,2))); %null correlation
end
disp('done!');
toc

%Ground truth correlations
A2 = A-mean(A,2); B2 = B-mean(B,2);
CCgt = sum(A2.*B2,2)./(sqrt(sum(A2.^2,2)).*sqrt(sum(B2.^2,2)));

%Z-score the correlation distributions
S = std(CCss(:));
CC_z = (CC-median(CC))./S;
CCss_z = (CCss-median(CCss,2))./S;
CCgt_z = (CCgt-median(CCgt))./S;

%compute histograms
edges = -10:0.2:10;
centers = (edges(1:end-1) + edges(2:end))/2;
hCC = histcounts(CC_z,edges);
hCCss = histcounts(CCss_z, edges); hCCss = (hCCss./sum(hCCss)).*sum(hCC);
hCCgt = histcounts(CCgt_z, edges); hCCgt = (hCCgt./sum(hCCgt)).*sum(hCC);

figure,
plot(centers, hCC, 'linewidth', 3)
hold on
plot(centers, hCCss, 'linewidth', 3)
plot(centers, hCCgt, 'linewidth', 3)
xlabel('Z score'); ylabel('Counts')

legend({'measured correlations', 'computed null', 'ground truth'});

[~,KSp] = kstest2(CC_z(:), CCss_z(:));

end