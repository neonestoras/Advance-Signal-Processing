
function pdf(x_in)

    num_bins = 100;
    width=max(x_in)-min(x_in);
    bin_size=width/num_bins;
    [counts,center] = hist(x_in,num_bins); 
    total_counts=sum(counts);
    rel_freq=counts./(total_counts*bin_size);
    bar(center,rel_freq,'FaceColor','b','EdgeColor','k')

end