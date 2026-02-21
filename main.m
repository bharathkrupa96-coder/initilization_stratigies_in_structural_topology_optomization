% initialisation study

for volfrac = [0.3 0.32 0.35 0.37 0.4 0.5 0.6 0.8 0.82 0.85 0.87 ]
    filename = ['comparison_results' num2str(volfrac) '.mat'];
    if ~isfile(filename)

        multiresults = Bop88_multistart_final(volfrac);
        save(filename)
    else
        load( filename);
    end
end




