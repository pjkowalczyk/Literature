################################################################
#BAYESIAN DECISION TREE
#Justin Klekota & Fritz Roth
#Harvard Medical School
#Department of Biological Chemistry and Molecular Pharmacology 
#250 Longwood Avenue, SGMB-322 Boston, MA 02115 
#Contact: klekota@gmail.com, froth@hms.harvard.edu
#Copyright 2008 President and Fellows of Harvard University 
################################################################

##########################################################################################################################
#INPUT FILE REQUIREMENTS:
#FILE must be a tab-delimited text file
#FIRST ROW contains column headers; each subsequent row contains data for one compound
#Total number of rows is less than $MAXLIBRARYSIZE
#FIRST COLUMN must be Compound ID
#SECOND COLUMN must be integer indicating activity (usually a 1 or 0)
#THIRD COLUMN must be integer indicating inactivity (usually a 1 or 0, opposite of second column)
#FOURTH AND HIGHER COLUMNS contain 1 or 0 indicating presence or absence of the corresponding substructure (SMILES/SMARTS)
##########################################################################################################################

##########################################################################################################################
#RUNTIME REQUIREMENTS:
#File::Sort Perl module required (CPAN.org)
#Execute Command: perl BayesianDecisionTree.pl YourDataFileName.txt
##########################################################################################################################

use strict;
use Sort qw(sort_file);

my $MAXLIBRARYSIZE=100000;

my $infilename=shift(@ARGV);
my $dl="\t";
my %Log=();
for(my $i=1; $i<$MAXLIBRARYSIZE; ++$i) {$Log{$i}=log($i);}

open(INFILE, $infilename) || die "Cannot open your infile $infilename\n";

my $headerString=<INFILE>;
chomp($headerString);
my @HeaderArray=split($dl, $headerString);

close STDOUT;
open(STDOUT, ">Tree.All.txt");
print STDOUT "This program assumes your compound ID and activity annotations are listed in the first ten columns\n";

my $max=10; if($max>scalar(@HeaderArray)) {$max=scalar(@HeaderArray);}
for(my $i=0; $i<$max; ++$i)
{
  print STDOUT $i.': '.$HeaderArray[$i]."\n";
}

print STDOUT "Which column number indicates the compound ID?\n";
my $numID=0; #<STDIN>; chomp($numID);
print STDOUT "Which column number indicates the number of active annotations?\n";
my $numActive=1; #<STDIN>; chomp($numActive);
print STDOUT "Which column number indicates the number of inactive annotations (this should NOT count untested)?\n";
my $numInactive=2; #<STDIN>; chomp($numInactive);
print STDOUT "Which column number indicates the first descriptor (substructure)?\n";
my $numBegin=3; #<STDIN>; chomp($numBegin);


#SORT INFILE BY COMPOUND ID
close INFILE;
File::Sort::sort_file({t => $dl, k => ($numID+1)."n", o => "temp.txt", I => $infilename});
open(INFILE, "temp.txt") || die "Cannot open temp.txt\n";
my $foo=<INFILE>; chomp($foo);
if($foo ne $headerString) {die "Header in sorted file did not stay at the top of the file.\n"
					       ."Make sure compound IDs are positive numbers.\n";}

#LOAD THE COMPOUND ACTIVITY/INACTIVITY DATA AND HASH THE COMPOUND IDs CONTAINING EACH SUBSTRUCTURE
my %Substructures=(); #stores arrays of Compound IDs containing each substructure tested
my %CompoundData=(); #contains the number of active and inactive annotations for each compound

#Decision Tree Objects
my %Node=(); #contains the path to each node of a given depth
my %CompoundList=(); #contains the compound IDs at a current leaf node
my %BIC=(); #contains the BIC score at each node;
my @Leaves=(); #contains the list of leaf nodes
my %NodeActivity=();
my %Significance=();
push(@{$Node{0}},""); #there is only one node of depth 0, i.e the root, which has an empty path ""

#Compound IDs should remain in order in all objects for the entire program
for(my $i=$numBegin; $i<scalar(@HeaderArray); ++$i)
{
    push(@{$Substructures{$HeaderArray[$i]}}, () );
}
while(my $tempString=<INFILE>)
{
  chomp($tempString);
  my @Temp=split($dl, $tempString);
  push(@{$CompoundData{$Temp[$numID]}},($Temp[$numActive],$Temp[$numInactive]));
  push(@{$CompoundList{""}}, $Temp[$numID]);
  if(scalar(@Temp)!=scalar(@HeaderArray)) {die "Header length mismatches Data Row length\n";}
  for(my $i=$numBegin; $i<scalar(@Temp); ++$i)
  {
    if($Temp[$i]) {push(@{$Substructures{$HeaderArray[$i]}}, $Temp[$numID]);}
  }
}
close INFILE;
unlink "temp.txt";

#Calculate BIC score at the root
my $active=0;
my $inactive=0;
my $N=0;
foreach my $tempCompound (@{$CompoundList{""}})
{
  $active+=$CompoundData{$tempCompound}[0];
  $inactive+=$CompoundData{$tempCompound}[1];
}
$N=$active+$inactive;
$BIC{""}=0;

my $rootP=$active/($active+$inactive);
if($active){$BIC{""}+=$active*log($active/$N);}
if($inactive){$BIC{""}+=$inactive*log($inactive/$N);}
$BIC{""}-=0.5*1*log($N); #N-1 parameters which is 2-1=1 in this case

print STDOUT "Sorted and loaded data successfully.\n\n";

#DECISION TREE BUILDING USING BIC SCORE
my $allFail=0;
my $currentDepth=0;
my $NodeCount=0;

while(!$allFail)
{
  print STDOUT "The current depth is ".$currentDepth."\n";
  my @CurrentDepthPaths=();
  push(@CurrentDepthPaths, @{$Node{$currentDepth}});
  $allFail=1; #Will be set back to 0 if split occurs so tree-building can continue
  for(my $zzz=0; $zzz<scalar(@CurrentDepthPaths); ++$zzz) #VISIT EACH NODE AT THIS DEPTH
  {
    close STDOUT;
    open(STDOUT, ">>Tree.All.txt");
    my $tempPath=$CurrentDepthPaths[$zzz];
    my $maxBIC=$BIC{$tempPath};
    my $maxSubIndex=-1;
    my $maxBICPlus=0;
    my $maxBICMinus=0;
    print STDOUT "The current depth is ".$currentDepth." testing node".$dl.$tempPath."\n";
    for(my $i=$numBegin; $i<scalar(@HeaderArray); ++$i) #RECALL THAT HEADER ARRAY CONTAINS SUBSTRUCTURE SMILES
    { 
      my $CC=0;
      my @PresentCompounds=();
      my @AbsentCompounds=();

      #Split Compounds at the Current Node into Two Groups With Given Substructure Present or Absent
      my $maxCC=scalar(@{$Substructures{$HeaderArray[$i]}});
      foreach my $tempCompound (@{$CompoundList{$tempPath}})
      {
	  while( (($Substructures{$HeaderArray[$i]}[$CC]) < $tempCompound) && ($CC<$maxCC) ) {++$CC;}
	  if($tempCompound == ($Substructures{$HeaderArray[$i]}[$CC]) ) {push(@PresentCompounds, $tempCompound);}
	  else {push(@AbsentCompounds, $tempCompound);}
      }

      my $tempActivePresent=0;
      my $tempInactivePresent=0;
      my $tempActiveAbsent=0;
      my $tempInactiveAbsent=0;
      my $tempBIC=0;      
      my $BICPlus=0;
      my $BICMinus=0;

      foreach my $tempCompound (@PresentCompounds)
      {
        $tempActivePresent+=$CompoundData{$tempCompound}[0];
        $tempInactivePresent+=$CompoundData{$tempCompound}[1];
      }
      foreach my $tempCompound (@AbsentCompounds)
      {
        $tempActiveAbsent+=$CompoundData{$tempCompound}[0];
        $tempInactiveAbsent+=$CompoundData{$tempCompound}[1];
      }
      #CALCULATE YOUR PSEUDOCOUNTS: USING OLIVER'S 2003 Genome Research Paper
      my $tempActives=$tempActivePresent+$tempActiveAbsent;
      my $tempInactives=$tempInactivePresent+$tempInactiveAbsent;
      my $adjustableParameter=1; #You can choose any number of pseudocounts, we choose 1
      my $pseudoCountActive=$adjustableParameter*$tempActives/($tempActives+$tempInactives);
      my $pseudoCountInactive=$adjustableParameter*$tempInactives/($tempActives+$tempInactives);
      $tempActivePresent+=$pseudoCountActive;
      $tempInactivePresent+=$pseudoCountInactive;
      $tempActiveAbsent+=$pseudoCountActive;
      $tempInactiveAbsent+=$pseudoCountInactive;

      my $NodeNPresent=$tempActivePresent+$tempInactivePresent;
      if($tempActivePresent>0)
      { $BICPlus+=$tempActivePresent*log($tempActivePresent/($tempActivePresent+$tempInactivePresent)); }
      if($tempInactivePresent>0)
      { $BICPlus+=$tempInactivePresent*log($tempInactivePresent/($tempActivePresent+$tempInactivePresent)); }
      $BICPlus-=0.5*1*log($N); #n-1 parameters which is 2-1=1 in this case

      my $NodeNAbsent=$tempActiveAbsent+$tempInactiveAbsent;
      if($tempActiveAbsent>0)
      { $BICMinus+=$tempActiveAbsent*log($tempActiveAbsent/($tempActiveAbsent+$tempInactiveAbsent)); }
      if($tempInactiveAbsent>0)
      { $BICMinus+=$tempInactiveAbsent*log($tempInactiveAbsent/($tempActiveAbsent+$tempInactiveAbsent)); }
      $BICMinus-=0.5*1*log($N); #n-1 parameters which is 2-1=1 in this case

      $tempBIC=$BICPlus+$BICMinus;
      if($tempBIC>$maxBIC)
      {
	$maxBIC=$tempBIC; $maxSubIndex=$i; $maxBICPlus=$BICPlus; $maxBICMinus=$BICMinus;
	my $UP="DOWN";
	if(($tempActivePresent/($tempInactivePresent+$tempActivePresent))>($tempActives/($tempActives+$tempInactives))) {$UP="UP";}
	@{$NodeActivity{$tempPath}}=($tempActives/($tempActives+$tempInactives),$tempActives,($tempActives+$tempInactives),
				     $tempActivePresent,$tempInactivePresent,$tempActiveAbsent,$tempInactiveAbsent,
				     $HeaderArray[$i], $UP, $NodeCount);
	++$NodeCount;
      }
      elsif($maxSubIndex<0)
      {
	@{$NodeActivity{$tempPath}}=($tempActives/($tempActives+$tempInactives),$tempActives,($tempActives+$tempInactives),
				     -1,-1,-1,-1,-1,-1, $NodeCount);
	++$NodeCount;
      }
    }
    if($maxSubIndex>=0)
    {
      $allFail=0; #indicates a change has been made so that the tree building will continue
      my $delta=$maxBIC-$BIC{$tempPath};
      #UPDATE DATA OBJECTS FOR THIS NEW SPLIT
      my $CC=0;
      my @PresentCompounds=();
      my @AbsentCompounds=();
      
      #Split Compounds at the Current Node into Two Groups With Given Substructure Present or Absent
      my $maxCC=scalar(@{$Substructures{$HeaderArray[$maxSubIndex]}});
      foreach my $tempCompound (@{$CompoundList{$tempPath}})
      {
	  while( (($Substructures{$HeaderArray[$maxSubIndex]}[$CC]) < $tempCompound) && ($CC<$maxCC) ) {++$CC;}
	  if($tempCompound == ($Substructures{$HeaderArray[$maxSubIndex]}[$CC]) ) {push(@PresentCompounds, $tempCompound);}
	  else {push(@AbsentCompounds, $tempCompound);}
      }

      @{$CompoundList{$tempPath}}=();
      my $plusPath=$tempPath.'+'.$dl.$HeaderArray[$maxSubIndex].$dl;
      my $minusPath=$tempPath.'-'.$dl.$HeaderArray[$maxSubIndex].$dl;
      push(@{$CompoundList{$plusPath}}, @PresentCompounds);
      push(@{$CompoundList{$minusPath}}, @AbsentCompounds);
      $BIC{$plusPath}=$maxBICPlus;
      $BIC{$minusPath}=$maxBICMinus;
      push(@{$Node{$currentDepth+1}}, ($plusPath, $minusPath));

      #####
      #MEASURE THE STATISTICAL SIGNIFICANCE OF THIS SPLIT HERE
      #($tempActives/($tempActives+$tempInactives),$tempActives,($tempActives+$tempInactives),
      #$tempActivePresent,$tempInactivePresent,$tempActiveAbsent,$tempInactiveAbsent);
      #
      my $Hit=$NodeActivity{$tempPath}[1];
      my $NonHit=$NodeActivity{$tempPath}[2]-$NodeActivity{$tempPath}[1];
      my $ClusterSizeTested=$NodeActivity{$tempPath}[3]+$NodeActivity{$tempPath}[4]-1;
      my $ClusterHitTemp=$NodeActivity{$tempPath}[3]-$NodeActivity{$tempPath}[0];
      
      my $hitMax=$ClusterSizeTested;
      if($ClusterSizeTested>$Hit) {$hitMax=$Hit;}

      #HYPERGEOMETRIC PROBABILITY H(I,C,H) (ACTIVITY)
      my $HICHTested=0;
      for(my $j=$ClusterHitTemp; $j<=$hitMax; ++$j)
      {
	  my $htemp=0;
	  for(my $k=($Hit-$j+1); $k<=$Hit; ++$k) {$htemp+=$Log{$k};}
	  for(my $k=1; $k<=$j; ++$k) {$htemp-=$Log{$k};}
	  for(my $k=($NonHit-$ClusterSizeTested+$j+1); $k<=($NonHit); ++$k) {$htemp+=$Log{$k};}
	  for(my $k=1; $k<=($ClusterSizeTested-$j); ++$k) {$htemp-=$Log{$k};}
	  for(my $k=($Hit+$NonHit-$ClusterSizeTested+1); $k<=$Hit+$NonHit; ++$k) {$htemp-=$Log{$k};}
	  for(my $k=1; $k<=$ClusterSizeTested; ++$k) {$htemp+=$Log{$k};}
	  $HICHTested+=exp($htemp);
      }
      my $hitMin=0;
      if($NonHit<$ClusterSizeTested) {$hitMin=$ClusterSizeTested-$NonHit;}

      #HYPERGEOMETRIC PROBABILITY H(I,C,N) (EXCLUSION)
      my $HICNTested=0;
      for(my $j=$hitMin; $j<=$ClusterHitTemp; ++$j)
      {
	  my $htemp=0;
	  for(my $k=($Hit-$j+1); $k<=$Hit; ++$k) {$htemp+=$Log{$k};}
	  for(my $k=1; $k<=$j; ++$k) {$htemp-=$Log{$k};}
	  for(my $k=($NonHit-$ClusterSizeTested+$j+1); $k<=($NonHit); ++$k) {$htemp+=$Log{$k};}
	  for(my $k=1; $k<=($ClusterSizeTested-$j); ++$k) {$htemp-=$Log{$k};}
	  for(my $k=($Hit+$NonHit-$ClusterSizeTested+1); $k<=$Hit+$NonHit; ++$k) {$htemp-=$Log{$k};}
	  for(my $k=1; $k<=$ClusterSizeTested; ++$k) {$htemp+=$Log{$k};}
	  $HICNTested+=exp($htemp);
      }
      if($HICNTested<$HICHTested) {$HICHTested=$HICNTested;}
      $Significance{$tempPath}=$HICHTested;
      #END SIGNIFICANCE
      #####

      print STDOUT "BIC score will increase ".$delta."\n";
      print STDOUT "The significance of this split is ".$Significance{$tempPath}." All $Hit of ".($NonHit+$Hit)." vs. Node $ClusterHitTemp of $ClusterSizeTested"."\n";
      print STDOUT "SPLIT! ".$NodeActivity{$tempPath}[0].': '.$NodeActivity{$tempPath}[8].', '.$HeaderArray[$maxSubIndex]."\n";
    }
    else
    {
      push(@Leaves, $tempPath);
      print STDOUT "LEAF! ".$NodeActivity{$tempPath}[0].': '.$tempPath."\n";
    }

  }
  ++$currentDepth;
}

close STDOUT;
open(STDOUT, ">>Tree.All.txt");

print STDOUT "\n\nLEAVES IN FINAL TREE!!!!!\n\n";
print STDOUT "The root node is ".$rootP." active\n\n";

for(my $i=0; $i<scalar(@Leaves); ++$i)
{
 my $tempLeaf=$Leaves[$i];
 print STDOUT $i.$dl;
 if($NodeActivity{$tempLeaf}[0]>$rootP) {print STDOUT "UP".$dl;} else {print STDOUT "DOWN".$dl;}
 print STDOUT $NodeActivity{$tempLeaf}[0].$dl.$NodeActivity{$tempLeaf}[1].$dl."of".$dl.$NodeActivity{$tempLeaf}[2].$dl.$tempLeaf."\n\n";
}

#@{$NodeActivity{$tempPath}}=($tempActives/($tempActives+$tempInactives),$tempActives,($tempActives+$tempInactives),
#                             $tempActivePresent,$tempInactivePresent,$tempActiveAbsent,$tempInactiveAbsent,
#			      $HeaderArray[$i], $UP, $NodeCount);

open(CL, ">CompoundLeaves.txt");
for(my $i=0; $i<scalar(@Leaves); ++$i)
{
  my $tempLeaf=$Leaves[$i];
  my $numCompounds=scalar(@{$CompoundList{$tempLeaf}});
  for(my $j=0; $j<$numCompounds; ++$j)
  {
    print CL $CompoundList{$tempLeaf}[$j].$dl.$i.$dl.$NodeActivity{$tempLeaf}[0]."\n";
  }
}
close CL;
