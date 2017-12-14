function AnnotatePolIIPeaks(transcripAnnotationsGFF, ...
                            peaksFilePrefix, ...
                            snapToAnnotationDist)
% Reads in the information about peaks in the PolII ChIP data, and
% tries to determine the direction of transcription for each PolII
% peak, using either transcript annotations from Ensembl or
% (failing that) the skewness of the distribution of fragment
% overlap counts in the peak.
% 
% Typical values would be:
%  transcripAnnotationsGFF = '../GenomeAnnotationsM37.gff';
%

peaksBedFileName = [peaksFilePrefix '.bed'];
peaksBarFileName = [peaksFilePrefix '.bar'];
startSitesBarFileName = [peaksFilePrefix '_StartSites.bar'];

[seqNames, ...
 sourceNames, ...
 featureNames, ...
 startCoords, ...
 endCoords, ...
 scores, ...
 strands] = textread(transcripAnnotationsGFF, ...
                        '%s %s %s %d %d %s %c %*[^\n]', ...
                        'delimiter', '\t', ...
                        'commentstyle', 'shell');
% GFF files are 1-based coordinate numbering; map to zero-based
startCoords = startCoords - 1;
endCoords = endCoords - 1;

chroms = unique(seqNames);
C = length(chroms);
seqNamesChrom = cell(1,C);
for c=1:C
  seqNamesChrom{c} = ['chr' chroms{c}];
end
chromAnnotInds = cell(1,C);
chromAnnotStartCoords = cell(1,C);
chromAnnotEndCoords = cell(1,C);
chromAnnotStrands = cell(1,C);
for c=1:C
  chromInds = find(strcmp(seqNames, chroms{c}));
  chromAnnotInds{c} = chromInds;
  chromAnnotStartCoords{c} = startCoords(chromInds);
  chromAnnotEndCoords{c} = endCoords(chromInds);
  chromAnnotStrands{c} = strands(chromInds);
end

[peakSeqNames, ...
 peakStartCoords, ...
 peakEndCoords, ...
 peakFeatNames, ...
 peakScores] = textread(peaksBedFileName, ...
                           '%s %d %d %s %f %*[^\n]', ...
                           'delimiter','\t', ...
                           'headerlines', 1);
% fix the "+1" end coordinate numbering in BED files
peakEndCoords = peakEndCoords - 1;

[barSeqNames, barData]=readbar(peaksBarFileName);

P = length(peakSeqNames);

outBarChroms = cell(P,1);
outBarCoords = zeros(P,1);
outBarData = zeros(P,1);

for p=1:P
  if mod(p, 1000) == 0
    fprintf(1, 'analyzed %d out of %d peaks\n', p, P);
  end
  
  % go through each peak
  peakChrom = peakSeqNames{p};
  peakStartCoord = peakStartCoords(p);
  peakEndCoord = peakEndCoords(p);
  
  barChromInd = find(strcmp(barSeqNames, peakChrom));
  if ~isempty(barChromInd)
    barChromData = barData{barChromInd};
    barChromCoords = barChromData(:,1);
    barDataInds = find(barChromCoords >= peakStartCoord & ...
                       barChromCoords <= peakEndCoord);
    barDataValues = barChromData(barDataInds,2);
    barDataCoords = barChromData(barDataInds,1);
    [peakMaxValue, maxValueInd] = max(barDataValues);
    peakCoord = barDataCoords(maxValueInd);
  else
    error('unable to get sample data for this peak');
  end  
  
  % which chromosome is the peak on?
  chromInd = find(strcmp(seqNamesChrom, peakChrom));
  if ~isempty(chromInd)
    % yes, we have Ensembl transcript annotations for this chromosome
    
    peakStrand = '';
    
    if ~isempty(chromInd)
      annotInds = chromAnnotInds{chromInd};
      annotStartCoords = chromAnnotStartCoords{chromInd};
      annotEndCoords = chromAnnotEndCoords{chromInd};
      annotStrands = chromAnnotStrands{chromInd};
      A = length(annotStartCoords);
      maxOverlap = 0;
      maxOverlapAnnotInd = [];
      overlapAnnotInd = [];
      for a=1:A
        fracOverlap = 0.5*(IntervalOverlap(peakStartCoord, ...
                                           peakEndCoord, ...
                                           annotStartCoords(a), ...
                                           annotEndCoords(a)) + ...
                           IntervalOverlap(annotStartCoords(a), ...
                                           annotEndCoords(a), ...
                                           peakStartCoord, ...
                                           peakEndCoord));
        
        if fracOverlap > maxOverlap
          
          % use the transcript annotation that maximally overlaps the PolII
          % peak, to determine the strand assignment of the PolII peak
          maxOverlap = fracOverlap;
          maxOverlapAnnotInd = annotInds(a);
          peakStrand = annotStrands(a);
          overlapAnnotInd = a;
        else
          if fracOverlap >= 0.4 && maxOverlap >= 0.4
            if annotStrands(a) ~= peakStrand
              % conflicting strand information! 
              peakStrand = '';
              overlapAnnotInd = [];
              break;
            end
          end
        end
      end
      if ~isempty(overlapAnnotInd) 
        if maxOverlap > 0.1
          % this peak region substantially overlaps an Ensembl
          % transcript annotation; use the 
          if annotStrands(overlapAnnotInd) == '+'
            annotPeakCoord = annotStartCoords(overlapAnnotInd);
          else
            annotPeakCoord = annotEndCoords(overlapAnnotInd);
          end
          if abs(peakCoord - annotPeakCoord) < snapToAnnotationDist
            peakCoord = annotPeakCoord;
          else
            % the start site for the transcript annotation is more
            % than 1 kb from the estimated peak; treat the
            % estimated peak location as definitive
          end
        else
          % this peak region only overlaps a transcript by a small
          % amount; it is more likely a novel transcription start
          % site or false positive, so don't use the transcript
          % annotation, and instead use the estimated 
        end
      end
    end      
  end

  outBarChroms{p} = peakChrom;
  outBarCoords(p) = peakCoord;
  outBarData(p) = peakMaxValue;

end

writebar(startSitesBarFileName, outBarChroms, outBarCoords, outBarData);
