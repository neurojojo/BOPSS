% Place all files in the folder below, then either specify which files contain
% the negative control data by inputting an array at line 6 below (or add 'neg'
% to the filename of the negative controls)

results = ResultsClass('C:\ying_12_18\ex556 D2\');
results.analyze % 
results.makeFigure
results.saveImages('show') % This will show the locations of puncta in the unoptimized search
results.saveImages('decoded','show') % This will show the locations of puncta in the optimized search

results.saveResults('file.xls')