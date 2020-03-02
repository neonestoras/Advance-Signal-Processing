
function pdf(x_in)

    num_bins = 100;
    bin_size=(max(x_in)-min(x_in))/num_bins; %samples in every bin
    [counts,centers] = hist(x_in,num_bins); 
    
    relative_freq=counts./(sum(counts)*bin_size); % convert counts to frequencY
    bar(centers,relative_freq,'FaceColor','b','EdgeColor','k')

end