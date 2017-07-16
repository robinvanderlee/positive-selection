#!/usr/bin/perl

package Guidance; #don't forget: a package must end with a return value (1; in the end)!!!!!

#use strict;
use FileHandle;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::Tools::CodonTable;
use File::Copy;
use FindBin qw($Bin);

use lib "$Bin/../bioSequence_scripts_and_constants/";
#use lib "/bioseq/bioSequence_scripts_and_constants/";



use MSA_parser;
use GENERAL_CONSTANTS;

use constant NEWIC2MAFFT =>              "$Bin/exec/newick2mafft.rb";
use constant MSA_SET_SCORE => 		     "$Bin/exec/msa_set_score";
use constant HOT_PROGRAM => 		     "$Bin/exec/HoT/COS.pl";
my $newick2mafft="$Bin/exec/newick2mafft.rb";
#my $MSA_Score_CSS="http://guidance.tau.ac.il/MSA_Colored.NEW.css";
my $MSA_Score_CSS="http://guidance.tau.ac.il/MSA_Colored.NEW.EM.css";

#my $phylonet_prog = "$Bin/exec/phylonet_v1_7/phylonet_v1_7.jar";
my $isEqualTopologyProg = "$Bin/../../programs/isEqualTree/isEqualTree";

sub CreateHTML_Graph
# Create an HTML BARs graph
# GET: 1. CSV FILE (The X var is the first Col and the Y is the second, X must be sorted)
#      2. HTML OUTPUT
#      3. X Lable (OPTIONAL)
############################################################################################
{
	my $CSV_File=shift; # The Data Must be in the correct X order
	my $Out=shift;
	my $X_LABLE=shift;
	open (OUT,">>$Out") || return ("Guidance::CreateHTML_Graph: Can't open Output '$Out' $!");
	open (DATA,$CSV_File) || return ("Guidance::CreateHTML_Graph: Can't open data file '$CSV_File' $!");
	# GRAPH PROPERTIES
	my $graphHeight = 200; # target height of graph
	my $BarWidth = 15; # width of bars
	my $maxResult = 1;
	my $scale = 1;
	my $last_x = 0;
	# Set the scale
	$scale = $graphHeight / $maxResult;
	my $line=<DATA>; # Header
	
	# SPACER BEFORE THE GRAPH
	print OUT "<tr><td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr>\n";
	print OUT "<tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr>\n";
	print OUT "<tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr>\n";
	
	print OUT "<tr>\n";
	print OUT "<td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\" style = 'text-align: right'>$X_LABLE<br>SCORE</td>\n";
	while ($line=<DATA>)
	{
		chomp ($line);
		my @data=split(",",$line);
#		print "$data[0]\t$last_x\n";
		if (($data[1] ne "NaN") and (($data[0]-$last_x)==1))
		{
			print OUT "<td valign = bottom style = 'border-bottom: 1px solid black";
			if ($data[0]==1)#The First Point, plot also the Y bar 
			{
				print OUT ";border-left: 1px solid black";
			}
			print OUT ";'><img title=\"$data[0]:$data[1]\" src=\"http://guidance.tau.ac.il/blue.gif\" width=\"$BarWidth\" height=\"".$data[1]*$scale."\" border=\"1\"></td>\n";
			$last_x=$data[0];
		}
		elsif ($data[1] ne "NaN") 
		{
			while ($data[0]-$last_x!=1)
			{
				print OUT "<td valign = bottom style = 'border-bottom: 1px solid black;'><img src=\"http://guidance.tau.ac.il/blue.gif\" width=\"$BarWidth\" height=\"0\" border=\"0\"></td>\n";
				$last_x++;
			}
			print OUT "<td valign = bottom style = 'border-bottom: 1px solid black;'><img title=\"$data[0]:$data[1]\" src=\"http://guidance.tau.ac.il/blue.gif\" width=\"$BarWidth\" height=\"".$data[1]*$scale."\" border=\"1\"></td>\n";
			$last_x=$data[0];
		}
		elsif ($data[1] eq "NaN")
		{
			print OUT "<td valign = bottom style = 'border-bottom: 1px solid black;'><img title=\"$data[0]:$data[1]\" src=\"http://guidance.tau.ac.il/blue.gif\" width=\"$BarWidth\" height=\"".$data[1]*$scale."\" border=\"1\"></td>\n";
			$last_x=$data[0];
		}
	}
	print OUT "</tr>\n";
	print OUT "</table>\n";
	print OUT "<b><P align=\"center\">Column</p></b>\n";
	close (OUT);
	close (DATA);
}

sub printColoredAlignment {
#############################################################################################################################################
#@ARGV == 3 or die "USAGE: $0 IN_MSA_FILE OUT_HTML_FILE SCORES_FILE
#SCORES_FILE - Each line should contain three values: column number, seq number, and a score between 0 and 1 (separated by white spaces)\n";
#############################################################################################################################################
	my $inMsaFile=shift;
	my $outHtmlFile=shift;
	my $scoresFile=shift;
	my $codesFile=shift; # OPTIONAL
	# Read scores
	open SCORES, $scoresFile or return ("Can't open $scoresFile: $!");
	my %scores;
	my %Code_Names;
	foreach (<SCORES>) {
		next if (/^#/);
		s/^\s+//;
		my ($col, $seq, $score) = split;
		$scores{$seq}[$col] = $score;
	}
	# Read Codes
	if ($codesFile ne "")
	{
		open (CODES,$codesFile) or return ("Guidance::printColoredAlignment Can't open the Codes file: '$codesFile' $!");
		while (my $line=<CODES>)
		{
			chomp $line;
			my ($Seq_name,$Code)=split("\t",$line);
			$Code_Names{$Code}=$Seq_name;
		}
		close (CODES);
	}
	# Read MSA
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $inMsaFile) or die "Can't open $inMsaFile: $!";
	my $aln = $in->next_aln;
	$aln->verbose(1); 
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat(); 
#	print "DEPTH:",$aln->num_sequences;<STDIN>;
	@ans=MSA_parser::check_msa_licit_and_size($inMsaFile,"fasta","no");
		if ($ans[0] eq "OK"){$MSA_Depth=$ans[1];}
		else {return "printColoredAlignment: ".join (" ",@ans);}
	# Print HTML start
	# Code from Conseq colored MSA: ~/pupkoSVN/trunk/www/conseq/runCalc_Conseq.pl line 985
	my %msaColors = ();
	my %msaPrintColors = ();
	my $lineCounter;
	my @line;
	my $key;

	my $fontSize=2;
	my $sequenceLengthForDisplay=400000;
	my @msaRightOrder=0;
	my $msaRightOrderCounter=0;
	my $tdWidth = 5;
	
	my @colorstep = (); #color steps
	$colorstep[0] = "#10C8D1"; #Not confident
	$colorstep[1] = "#8CFFFF";
	$colorstep[2] = "#D7FFFF";
	$colorstep[3] = "#EAFFFF";
	$colorstep[4] = "#FFFFFF"; #average
	$colorstep[5] = "#FCEDF4";
	$colorstep[6] = "#FAC9DE";
	$colorstep[7] = "#F07DAB";
	$colorstep[8] = "#A02560"; #Most confident
	$colorstep[9] = "#A02560"; #Most confident (the score is exactly 1)
	open MSACOLOREDHTML, ">$outHtmlFile" or die "Can't open $outHtmlFile: $!";
	print MSACOLOREDHTML "<html>\n<head>\n</head>\n<body>\n\n";
	print MSACOLOREDHTML "<H1 align=center><u>MSA color-coded by GUIDANCE scores</u></H1>\n\n";
	print MSACOLOREDHTML "<table border=0  CELLSPACING=1  CELLPADDING=0 >\n";


	# Print colored HTML
	
	# counts how many times we print the whole section (relevants to sequences longer than the sequenceLengthForDisplay) 
	for(my $blockStart=1; $blockStart<$aln->length; $blockStart+=$sequenceLengthForDisplay) {
		my $blockEnd = $blockStart+$sequenceLengthForDisplay;
		$blockEnd = $aln->length if ($blockEnd > $aln->length);

		# Iterate over sequences and print up to sequenceLengthForDisplay residues
#		foreach my $seq ($aln->each_seq) { #HAIM COMMENT
#		my $depth = $aln->num_sequences;	    #HAIM ADD
		for(my $i=1;$i<=$MSA_Depth;$i++){		    #HAIM ADD
			my $seq = $aln->get_seq_by_pos($i); #HAIM ADD

			#		next if (   $seq->id < 42
			#				 || $seq->id > 73);
	
			# Print seq id
			print MSACOLOREDHTML "<tr>\n";
			print MSACOLOREDHTML "<td><b><font face='Courier New' color='black' size=$fontSize>", $seq->id, "</font></b></td>\n" if ($codesFile eq "");
			print MSACOLOREDHTML "<td><b><font face='Courier New' color='black' size=$fontSize>", substr($Code_Names{$seq->id},0,25), "</font></b></td>\n" if ($codesFile ne "");	
	
			# Print seq
			my @seq = split //, $seq->subseq($blockStart, $blockEnd);

			for(my $pos=0; $pos<@seq; $pos++) {
				my $res = $seq[$pos];
				if ($res eq '-') {
					print MSACOLOREDHTML "<td width=$tdWidth><b><font face='Courier New' color='black' size=$fontSize>$res</font></b></td>\n";
				} else {
					#print $seq->id,"\tscores{$seq->id}[$pos+1]:$scores{$seq->id}[$pos+1]\n";
					#my $color = $colorstep[ int(9 * $scores{$seq->id}[$pos+1]) ]; #HAIM COMMENT
					my $color = $colorstep[ int(9 * $scores{$i}[$pos+1]) ];
					if($color eq "#A02560"){
						print MSACOLOREDHTML "<td width=$tdWidth><b><font face='Courier New' color='white' size=$fontSize><span style='background: $color;'>$res</span></font></b></td>\n";
					}
					else {                
						print MSACOLOREDHTML "<td width=$tdWidth><b><font face='Courier New' color='black' 	size=$fontSize><span style='background: $color;'>$res</span></font></b></td>\n";
					}
				}
			}
			print MSACOLOREDHTML "</tr>\n\n";
		}
		
		print MSACOLOREDHTML "<tr>&nbsp</tr>\n\n";
	}
	print MSACOLOREDHTML "</table>";
		
	# print the color scale
	print  MSACOLOREDHTML "\n<br><b><u>Legend:</u><br><br>\nThe alignment confidence scale:</b><br>\n<table border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
	for (my $i=8 ; $i>=0 ; $i--){
		if ($i == 0){
			print  MSACOLOREDHTML "<font face='Courier New' color='white' size=$fontSize><span style='background: $colorstep[$i];'>&nbsp;", $i+1, "&nbsp;</span></font>";
		}
		else {
			print  MSACOLOREDHTML "<font face='Courier New' color='black' size=$fontSize><span style='background: $colorstep[$i];'>&nbsp;", $i+1, "&nbsp;</span></font>";
		}
	}
	print  MSACOLOREDHTML "</font></center>\n<center><table border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Confident</b></td>\n<td align=center><b><---></b></td>\n<td align=right><b>Uncertain</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
	print MSACOLOREDHTML "</body>\n<html>\n";
	close MSACOLOREDHTML;
	close (SCORES);
}

sub AssignColorsToAlignment{
	my $inMsaFile=shift;
	my $scoresFile=shift;
	my $codesFile=shift; # OPTIONAL

	open SCORES, $scoresFile or return ("Can't open $scoresFile: $!");
	my %scores;
	my %Code_Names;
	foreach (<SCORES>) {
		next if (/^\#/);
		s/^\s+//;
		my ($col, $seq, $score) = split;
		$scores{$seq}[$col] = $score;
	}
	# Read Codes
	if ($codesFile ne "")
	{
		open (CODES,$codesFile) or return ("Guidance::printColoredAlignment Can't open the Codes file: '$codesFile' $!");
		while (my $line=<CODES>)
		{
			chomp $line;
			my ($Seq_name,$Code)=split("\t",$line);
			$Code_Names{$Code}=$Seq_name;
		}
		close (CODES);
	}
	# Read MSA
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $inMsaFile) or die "Can't open $inMsaFile: $!";
	my $aln = $in->next_aln;
	$aln->verbose(1);
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat();
#	print "DEPTH:",$aln->num_sequences;<STDIN>;
	@ans=MSA_parser::check_msa_licit_and_size($inMsaFile,"fasta","no"); #HAIM ADD
	if ($ans[0] eq "OK"){$MSA_Depth=$ans[1];}    #HAIM ADD
	else {return "printColoredAlignment: ".join (" ",@ans);} #HAIM ADD
	# Print HTML start
	# Code from Conseq colored MSA: ~/pupkoSVN/trunk/www/conseq/runCalc_Conseq.pl line 985
	my %msaColors = ();
	my %msaPrintColors = ();
	my $lineCounter;
	my @line;
	my $key;
	
	my $fontSize=2;
	my $sequenceLengthForDisplay=400000;
	my @msaRightOrder=0;
	my $msaRightOrderCounter=0;
	my $tdWidth = 5;
	
	my @colorstep = (); #color steps
	$colorstep[0] = "Score1"; #Not confident
	$colorstep[1] = "Score2";
	$colorstep[2] = "Score3";
	$colorstep[3] = "Score4";
	$colorstep[4] = "Score5"; #average
	$colorstep[5] = "Score6";
	$colorstep[6] = "Score7";
	$colorstep[7] = "Score8";
	$colorstep[8] = "Score9"; #Most confident
	$colorstep[9] = "Score9"; #Most confident (the score is exactly 1)
	# get Align max_seq_length
	my $seq = $aln->get_seq_by_pos(1);
	my $Align_width = $seq->length();
	# STOPPED HERE
	

}

sub printColoredAlignment_With_CSS {
#############################################################################################################################################
#@ARGV == 3 or die "USAGE: $0 IN_MSA_FILE OUT_HTML_FILE SCORES_FILE
#SCORES_FILE - Each line should contain three values: column number, seq number, and a score between 0 and 1 (separated by white spaces)\n";
#############################################################################################################################################
	my $inMsaFile=shift;
	my $outHtmlFile=shift;
	my $scoresFile=shift;
	my $codesFile=shift; # OPTIONAL
#	my $COL_SCORES_FIGURE=shift; # OPTIONAL
#	my $MSA_Length=shift;
    # Parameters for the PLOT beneath alignment
	my $ColScoresCSV=shift;
	my $XLable=shift;
	my $Seq_Scores=shift;
	# Read scores
	open (SEQ_SCORES,$Seq_Scores) or return ("Can't open $Seq_Scores: $!");
	my %seq_scores=();
	foreach (<SEQ_SCORES>)
	{
		next if (/^#/);
		s/^\s+//;
        my ($seq, $score) = split;
        $seq_scores{$seq} = $score;
	}
		close (SEQ_SCORES);
	open SCORES, $scoresFile or return ("Can't open $scoresFile: $!");
	my %scores;
	my %Code_Names;
	foreach (<SCORES>) {
		next if (/^#/);
		s/^\s+//;
		my ($col, $seq, $score) = split;
		$scores{$seq}[$col] = $score;
	}
	# Read Codes
	if ($codesFile ne "")
	{
		open (CODES,$codesFile) or return ("Guidance::printColoredAlignment Can't open the Codes file: '$codesFile' $!");
		while (my $line=<CODES>)
		{
			chomp $line;
			my ($Seq_name,$Code)=split("\t",$line);
			$Code_Names{$Code}=$Seq_name;
		}
		close (CODES);
	}
	# Read MSA
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $inMsaFile) or die "Can't open $inMsaFile: $!";
	my $aln = $in->next_aln;
	$aln->verbose(1); #HAIM COMMNET
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat(); #HAIM COMMENT
#	print "DEPTH:",$aln->num_sequences;<STDIN>;
	@ans=MSA_parser::check_msa_licit_and_size($inMsaFile,"fasta","no"); #HAIM ADD
	if ($ans[0] eq "OK"){$MSA_Depth=$ans[1];}    #HAIM ADD
	else {return "printColoredAlignment: ".join (" ",@ans);} #HAIM ADD
	# Print HTML start
	# Code from Conseq colored MSA: ~/pupkoSVN/trunk/www/conseq/runCalc_Conseq.pl line 985
	my %msaColors = ();
	my %msaPrintColors = ();
	my $lineCounter;
	my @line;
	my $key;

	my $fontSize=2;
	my $sequenceLengthForDisplay=400000;
	my @msaRightOrder=0;
	my $msaRightOrderCounter=0;
	my $tdWidth = 5;
	
	my @colorstep = (); #color steps
	$colorstep[0] = "Score1"; #Not confident
	$colorstep[1] = "Score2";
	$colorstep[2] = "Score3";
	$colorstep[3] = "Score4";
	$colorstep[4] = "Score5"; #average
	$colorstep[5] = "Score6";
	$colorstep[6] = "Score7";
	$colorstep[7] = "Score8";
	$colorstep[8] = "Score9"; #Most confident
	$colorstep[9] = "Score9"; #Most confident (the score is exactly 1)
	
	my @colorstep_code = (); #color steps
	$colorstep_code[0] = "#10C8D1"; #Not confident
	$colorstep_code[1] = "#8CFFFF";
	$colorstep_code[2] = "#D7FFFF";
	$colorstep_code[3] = "#EAFFFF";
	$colorstep_code[4] = "#FFFFFF"; #average
	$colorstep_code[5] = "#FCEDF4";
	$colorstep_code[6] = "#FAC9DE";
	$colorstep_code[7] = "#F07DAB";
	$colorstep_code[8] = "#A02560"; #Most confident
#	my $plot_width=0;
#	if ($MSA_Length>1000)
#	{
#		   $plot_width=11.15*$MSA_Length;
#	   }
#		else
#	{
#		$plot_width=11.5*$MSA_Length;
#	   }
		open MSACOLOREDHTML, ">$outHtmlFile" or die "Can't open $outHtmlFile: $!";
		print MSACOLOREDHTML "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\"\n";
		print MSACOLOREDHTML "\"http://www.w3.org/TR/html4/strict.dtd\">\n";
		print MSACOLOREDHTML "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"\n";
		print MSACOLOREDHTML "\"http://www.w3.org/TR/html4/loose.dtd\">\n";
		print MSACOLOREDHTML "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Frameset//EN\"\n";
		print MSACOLOREDHTML "\"http://www.w3.org/TR/html4/frameset.dtd\">\n";
		print MSACOLOREDHTML "<head>\n";
		print MSACOLOREDHTML "<meta http-equiv=\"X-UA-Compatible\" content=\"IE=EmulateIE7\"/>\n";
#	print MSACOLOREDHTML "<html>\n";
#	print MSACOLOREDHTML "<head>\n";
	print MSACOLOREDHTML "<link rel=\"stylesheet\" type=\"text/css\" href=\"$MSA_Score_CSS\"/>\n";
#	print MSACOLOREDHTML "<style type=\"text/css\">\n";
#	print MSACOLOREDHTML "img.plot{\n";
#   print MSACOLOREDHTML "width: $plot_width"."px;\n";
##	print MSACOLOREDHTML "border-right: 2px dotted #4169e1;\n";
##	print MSACOLOREDHTML "border-bottom: 2px dotted #4169e1;\n";
#    print MSACOLOREDHTML "}\n";
#	print MSACOLOREDHTML "</style>\n";

	print MSACOLOREDHTML "</head>\n";
	print MSACOLOREDHTML "<H1 align=center><u>MSA color-coded by GUIDANCE scores</u></H1>\n\n";
	print MSACOLOREDHTML "<table>\n";
	# Print colored HTML
	
	# get Align max_seq_length
        my $seq = $aln->get_seq_by_pos(1);
        my $Align_width = $seq->length();
	
        # Print upper Scale
	print MSACOLOREDHTML "<tr>\n<td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</td><td>1</td>\n";
	my $i=2;
        while ($i<$Align_width)
#	for (my $i=2;$i<$Align_width;$i++)
	  {
	    if (($i%10)==0)
	      {
		my @digits=split("",$i);
		print MSACOLOREDHTML "\n";
		foreach my $digit (@digits)
		  {
		    print MSACOLOREDHTML "<td>$digit</td>";
		    $i++;
		  }
		
	      }
	    else
	      {
		print MSACOLOREDHTML "<td></td>";
		$i++;
	      }
	  }
	  print  MSACOLOREDHTML "\n</tr>\n";

	# counts how many times we print the whole section (relevants to sequences longer than the sequenceLengthForDisplay) 
	for(my $blockStart=1; $blockStart<$aln->length; $blockStart+=$sequenceLengthForDisplay) {
		my $blockEnd = $blockStart+$sequenceLengthForDisplay;
		$blockEnd = $aln->length if ($blockEnd > $aln->length);

		# Iterate over sequences and print up to sequenceLengthForDisplay residues
#		foreach my $seq ($aln->each_seq) { #HAIM COMMENT
#		my $depth = $aln->num_sequences;	    #HAIM ADD
		for(my $i=1;$i<=$MSA_Depth;$i++){		    #HAIM ADD
			my $seq = $aln->get_seq_by_pos($i); #HAIM ADD

	#		next if (   $seq->id < 42
	#				 || $seq->id > 73);
	
			# Print seq id
			print MSACOLOREDHTML "<tr>\n";
			# print SEQ COLOR
			my $Color_Class="";
			if ($seq_scores{$i} eq "nan") {$Color_Class="ScoreNaN";}
			else {$Color_Class = $colorstep[ int(9 * $seq_scores{$i})];}
            print MSACOLOREDHTML "<td class=\"$Color_Class\">&nbsp</td>\n";
			print MSACOLOREDHTML "<td class=\"Seq_Name\">", $seq->id, "</td>\n" if ($codesFile eq "");
			print MSACOLOREDHTML "<td class=\"Seq_Name\">", substr($Code_Names{$seq->id},0,25), "</td>\n" if ($codesFile ne "");	
				
			# Print seq
			my @seq = split //, $seq->subseq($blockStart, $blockEnd);

			for(my $pos=0; $pos<@seq; $pos++) {
	#		for(my $pos=757; $pos<875; $pos++) {
				my $res = $seq[$pos];
				if ($res eq '-') {
					print MSACOLOREDHTML "<td>$res</td>\n";
				} else {
					#print $seq->id,"\tscores{$seq->id}[$pos+1]:$scores{$seq->id}[$pos+1]\n";
					#my $color = $colorstep[ int(9 * $scores{$seq->id}[$pos+1]) ]; #HAIM COMMENT
				  	my $Color_Class="";
					if ($scores{$i}[$pos+1] eq "nan") {$Color_Class="ScoreNaN";}
					else {$Color_Class = $colorstep[ int(9 * $scores{$i}[$pos+1]) ]};
					print MSACOLOREDHTML "<td class=\"$Color_Class\">$res</td>\n";
				}
			}
			print MSACOLOREDHTML "</tr>\n\n";
		}

#		print MSACOLOREDHTML "<tr>&nbsp</tr>\n\n";
  	    print MSACOLOREDHTML "<tr></tr>\n\n";
	}
	

	
        
      # Print lower Scale
#	print MSACOLOREDHTML "<tr>\n<td class=\"Seq_Name\">&nbsp</td><td>1</td>\n";
	print MSACOLOREDHTML "<tr>\n<td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\"></td><td>1</td>\n";
	$i=2;
        while ($i<$Align_width)
#	for (my $i=2;$i<$Align_width;$i++)
	  {
	    if (($i%10)==0)
	      {
		my @digits=split("",$i);
		foreach my $digit (@digits)
		  {
		    print MSACOLOREDHTML "\n<td>$digit</td>";
		    $i++;
		  }
		
	      }
	    else
	      {
		print MSACOLOREDHTML "<td></td>";
		$i++;
	      }
	  }
        print  MSACOLOREDHTML "\n</tr>\n";
						
#		print MSACOLOREDHTML "<tr>&nbsp</tr>\n</table>\n";

        # Add the figure
		close (MSACOLOREDHTML);
		CreateHTML_Graph($ColScoresCSV,$outHtmlFile,$XLable);
		open MSACOLOREDHTML, ">>$outHtmlFile";
#        if ($COL_SCORES_FIGURE ne "")
#	  {
#            print MSACOLOREDHTML "<table>\n";
#	    print MSACOLOREDHTML "<tr>\n";
#	    print MSACOLOREDHTML "<td class=\"Seq_Name\">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</td>\n";
#	    print MSACOLOREDHTML "<td><img class=\"plot\" src=\"$COL_SCORES_FIGURE\"></td>\n";
#	    print MSACOLOREDHTML "</table>\n";
#	  }
						
	# print the color scale
	print  MSACOLOREDHTML "\n<br><b><u>Legend:</u><br><br>\nThe alignment confidence scale:</b><br>\n";#<table border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
		print MSACOLOREDHTML "<table style = 'table-layout: auto;margin-left: 0em;margin-right: 0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse: collapse;' border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
	for (my $i=8 ; $i>=0 ; $i--){
    
		if ($i == 0){
		
			print  MSACOLOREDHTML "<font face='Courier New' color='white' size=$fontSize><span style='background: $colorstep_code[$i];'>&nbsp;", $i+1, "&nbsp;</span></font>";
		}
		else {
		
			print  MSACOLOREDHTML "<font face='Courier New' color='black' size=$fontSize><span style='background: $colorstep_code[$i];'>&nbsp;", $i+1, "&nbsp;</span></font>";
		}
	}
#	print  MSACOLOREDHTML "</font></center>\n<center><table border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Confident</b></td>\n<td align=center><b><---></b></td>\n<td align=right><b>Uncertain</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
		print  MSACOLOREDHTML "</font></center>\n<center><table style = 'table-layout: auto;margin-left:0em;margin-right: 0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse: collapse;' border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Confident</b></td>\n<td align=center><b><---></b></td>\n<td align=right><b>Uncertain</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";

#    print MSACOLOREDHTML "<left><table border=0 cols=2 width=100>\n<tr>\n<td align=center class=\"ScoreNaN\">&nbsp;</td><td align=left>Insufficient Data</b></td>";
     print MSACOLOREDHTML "<left><table style = 'table-layout: auto;margin-left: 0em;margin-right:0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse:collapse;' border=0 cols=2 width=100>\n<tr>\n<td align=center class=\"ScoreNaN\">&nbsp;</td><td align=left>Insufficient Data</b></td>";

        
	print MSACOLOREDHTML "</body>\n</table>\n";
	close MSACOLOREDHTML;
}

sub root_BP_trees
####################################################################################################################
# Root all trees on BP dir
####################################################################################################################
{
	my $bsDir = shift;
	my $dataset = shift;
	my $orig_prog = shift;
	my $bp_repeats = shift;
	my $suffix=shift;
	if (!defined $suffix) {
    		$suffix = "";
	}
	$bsDir .= "/";

	for (my $countTrees=0;$countTrees<$bp_repeats;++$countTrees) {
		#print "**** $dataset tree_$countTrees\n";
		my $treeFile = $bsDir."tree_".$countTrees."/".$dataset.".".$orig_prog.".semphy.tree_".$countTrees.$suffix;
		if (-e $treeFile) {
			my $rootedTreeFile = $treeFile.".rooted";
#        		print "treeFile: $treeFile\nrootedTreeFile: $rootedTreeFile\n";
			rootTree($treeFile,$rootedTreeFile);
			open IN, "<$rootedTreeFile", or return ("can't open $rootedTreeFile");
			my $newick = "";
			foreach my $line (<IN>) {
				$newick.=$line;
			}
			close IN;
			
			#       if ($newick =~ m/:0[^\.]/ && $newick !~ m/:0\./) {
			if ($newick =~ m/:-/) {
				my $rootedTreeFile_withMinusLengths = $rootedTreeFile.".withMinusLengths";
				my $command = "mv $rootedTreeFile $rootedTreeFile_withMinusLengths";
				system ($command);
				$newick =~ s/:-/:/g;
				$command = "echo '$newick' > $rootedTreeFile";
				system ($command);
			}
			if ($newick =~ m/:\d+[^\.\d+]/) {
				my $rootedTreeFile_badBranchLength = $rootedTreeFile.".badBranchLength";
				my $command = "mv $rootedTreeFile $rootedTreeFile_badBranchLength";
				system ($command);
				$newick =~ s/:(\d+)([^\.\d+])/:$1\.0$2/g;
				# print "$newick\n";
				$command = "echo '$newick' > $rootedTreeFile";
				system ($command);
			}
		}else {
			print "File does not exist: $treeFile\n";
		}
	}
	return "ok";
}
sub convertBPTrees2MafftFormat{
####################################################################################################################
#convert the tree to mafft format
####################################################################################################################
#die " $0 <BPdir> <dataset> <original alignment program> <num_BP_repeats> <extra suffix - default ''> required as input\n" if (scalar(@ARGV) < 4);

	my $bsDir = shift;
	my $dataset = shift;
	my $orig_prog = shift;
	my $bp_repeats = shift;
	my $ruby_path=shift;
	my $suffix=shift;
	if (!defined $ruby_path)
	{
		$ruby_path="ruby";
	}
	if (!defined $suffix) {
    		$suffix = "";
	}
	
	$bsDir .= "/";

	for (my $countTrees=0;$countTrees<$bp_repeats;++$countTrees) {
    		#print "**** $dataset tree_$countTrees\n";
    		my $treeFile = $bsDir."tree_".$countTrees."/".$dataset.".".$orig_prog.".semphy.tree_".$countTrees.$suffix;
    		if (-e $treeFile) {
        		my $rootedTreeFile = $treeFile.".rooted";
        		my $mafftFormatTreeFile = $treeFile.".mafftFormat";
#        		print "treeFile: $treeFile\nrootedTreeFile: $rootedTreeFile\nmafftFormatTreeFile: $mafftFormatTreeFile\n";
		
        		rootTree($treeFile,$rootedTreeFile);
        		open IN, "<$rootedTreeFile", or return ("can't open $rootedTreeFile");
        		my $newick = "";
        		foreach my $line (<IN>) {
            			$newick.=$line;
        		}
        		close IN;

		#       if ($newick =~ m/:0[^\.]/ && $newick !~ m/:0\./) {
	        	if ($newick =~ m/:-/) {
            			my $rootedTreeFile_withMinusLengths = $rootedTreeFile.".withMinusLengths";
            			my $command = "mv $rootedTreeFile $rootedTreeFile_withMinusLengths";
            			system ($command);
            			$newick =~ s/:-/:/g;
            			$command = "echo '$newick' > $rootedTreeFile";
            			system ($command);
        		}
        		if ($newick =~ m/:\d+[^\.\d+]/) {
	            		my $rootedTreeFile_badBranchLength = $rootedTreeFile.".badBranchLength";
            			my $command = "mv $rootedTreeFile $rootedTreeFile_badBranchLength";
            			system ($command);
            			$newick =~ s/:(\d+)([^\.\d+])/:$1\.0$2/g;
				# print "$newick\n";
            			$command = "echo '$newick' > $rootedTreeFile";
            			system ($command);
        		}

        		system("$ruby_path $newick2mafft $rootedTreeFile > $mafftFormatTreeFile\n");
    		}else {
        		print "File does not exist: $treeFile\n";
    		}
	}
	return "ok";
}


sub runAlignBPtrees{
#die " $0 <alignment program (prank/mafft/clustal/muscle/pagan)> <dataset> <noBPdir> <original alignment program - default = same as current program> <num_BP_repeats - default 20> <extra suffix - default ''> required as input\n" if (scalar(@ARGV) < 3);
	my $aln_program = "";
	my $prog = shift;

	my $dataset = shift;
	my $noBPdir = shift;
	if ($noBPdir !~/(\/)$/){$noBPdir.="/";}
	my $aminoFile = shift;
	
	my $Seq_Type = shift;
	if ($Seq_Type eq "") {$Seq_Type="AminoAcids";}
	my $orig_prog=shift;
	if ($orig_prog eq "") {$orig_prog = $prog;}
	my $bp_repeats = shift;
	if ($bp_repeats eq "") {$bp_repeats = 20;}
	my $align_param=shift; # alignment additional parameters
	if (!defined $align_param) {$align_param = "";}
	my $suffix = shift;
	if (!defined $suffix) {$suffix = "";}
	my $update_file=shift; #OPTIONAL FOR SERVER TO UPDATE THE PROGRESS
	my $mafft_prog = shift;
	my $prank_prog = shift;
	my $clustalw_prog = shift;
	my $muscle_prog = shift;
	my $pagan_prog=shift;
	my $proc_num = shift; #AMIT

	my $PRANK_VERSION=0;
	if (uc($prog) eq "PRANK")  {
		if ($Seq_Type eq "AminoAcids") {$aln_program = $prank_prog;}
		elsif ($Seq_Type eq "Nucleotides") {$aln_program = $prank_prog;}
		elsif ($Seq_Type eq "Codons") {$aln_program = $prank_prog." -translate ";}
		# if PRANK find out it version
		my $prank_help=`$aln_program`;
		if ($prank_help=~/.*prunedata.*/){$PRANK_VERSION="121218";}
		elsif ($prank_help=~/.*showanc.*/){$PRANK_VERSION="120626";}
		#print "PRANK VERSION:>=$PRANK_VERSION\n";
	} elsif (uc($prog) eq "MAFFT"){
		if ($Seq_Type eq "AminoAcids") {$aln_program = $mafft_prog." --quiet --amino  ";}
		elsif ($Seq_Type eq "Nucleotides") {$aln_program = $mafft_prog." --quiet --nuc ";}
	} elsif (uc($prog) eq "MAFFT_LINSI"){
		if ($Seq_Type eq "AminoAcids") {$aln_program = GENERAL_CONSTANTS::MAFFT_LINSI_GUIDANCE." --amino  ";}
		elsif ($Seq_Type eq "Nucleotides") {$aln_program = GENERAL_CONSTANTS::MAFFT_LINSI_GUIDANCE." --nuc ";}
	} elsif (uc($prog) eq "CLUSTALW"){
		if ($Seq_Type eq "AminoAcids") {$aln_program = $clustalw_prog." -QUIET -TYPE=PROTEIN ";}
		elsif ($Seq_Type eq "Nucleotides") {$aln_program = $clustalw_prog." -QUIET -TYPE=DNA ";}
	} elsif (uc($prog) eq "MUSCLE"){
		if ($Seq_Type eq "AminoAcids") {$aln_program = $muscle_prog." -quiet -seqtype protein ";}
		elsif ($Seq_Type eq "Nucleotides") {$aln_program = $muscle_prog." -quiet -seqtype dna ";}
	} elsif (uc($prog) eq "PAGAN"){
		$aln_program=$pagan_prog." ";
	}
 	
	else {
    		return "please specify either mafft / mafft_linsi / prank / clustalw / muscle / pagan as an alignment program\n";
	}
	my $dir = $noBPdir."BP/";
#	my $aminoFile = $noBPdir.$dataset.".fas".$suffix;

	# Running parallel processes using fork, each will run an equal share of the BP alignments 
#	my $proc_num = 8;
	my $bp_per_proc = int($bp_repeats / $proc_num) + 1;
	my @children;  # array of pid of children
	for (my $proc=0; $proc<$proc_num; $proc++) {
		print "Running proc num $proc\n";

		my $pid = fork();
		if ($pid) {

			# parent
			push(@children, $pid);

		} elsif ($pid == 0) {

			# child
			for (my $tree_num=0; $tree_num<$bp_per_proc; $tree_num++) {
				my $countTrees = $proc * $bp_per_proc + $tree_num;
				print "proc num $proc\ttree num $tree_num --> global tree index $countTrees\n";
				last if ($countTrees == $bp_repeats);
				my $treeFile = $dir."tree_$countTrees/".$dataset.".".$orig_prog.".semphy.tree_".$countTrees.$suffix;
				if (-e $treeFile) {
					my $alnFile = $dir."tree_$countTrees/".$dataset.".".$orig_prog.".tree_$countTrees.".$prog.".aln";
					my $cmdFile = $dir."tree_$countTrees/".$dataset.".".$orig_prog.".tree_$countTrees.".$prog.".cmd";
					my $stdFile = $dir."tree_$countTrees/".$dataset.".".$orig_prog.".tree_$countTrees.".$prog.".std";
					#	if (-e $alnFile) {
					#	    print "skipping tree $countTrees because $alnFile already exists.\n";
					#	    next;
					#	}	
					my $cmd = "";
					if (uc($prog) eq "PRANK") {
						if ($PRANK_VERSION==0)
						{
							$cmd = "$aln_program $align_param -d=$aminoFile -o=$alnFile -t=$treeFile -once -quiet -noxml > $stdFile";# >& $stdFile"; # -noxml is not supported any more
						}
						else # version >=120626
						{
							$cmd = "$aln_program $align_param -d=$aminoFile -o=$alnFile -t=$treeFile -once -quiet > $stdFile"; # -noxml is not supported any more
						}
					} elsif ((uc($prog) eq "MAFFT") or (uc($prog) eq "MAFFT_LINSI")){
						my $mafftFormatTreeFile = $treeFile.".mafftFormat";
						unless (-e $mafftFormatTreeFile) {
							return ("ERROR: file does not exist: $mafftFormatTreeFile\n");
						}
						if ($align_param=~/addfragments/)
						{
							# first build the core MSA
							my $core_aln_param="";
							my @tmp=split (/\Q\-\-\E/,$align_param);
							my $tmp_size=@tmp;
							if ($tmp_size>=1) # command line type
							{
								for (my $i=0;$i<$tmp_size;$i++)
								{
									if ($tmp[$i]=~/addfragments/){delete $tmp[$i];}
									elsif ($tmp[$i]=~/multipair/){delete $tmp[$i];}
									elsif ($tmp[$i]=~/6merpair/){delete $tmp[$i];}
								}
								$core_aln_param=join ("\\\-\\\-",@tmp);
							}
							else # server type
							{
								my @tmp=split (/\Q--\E/,$align_param);
								my $tmp_size=@tmp;
								for (my $i=0;$i<$tmp_size;$i++)
								{
									if ($tmp[$i]=~/addfragments/){delete $tmp[$i];}
									elsif ($tmp[$i]=~/multipair/){delete $tmp[$i];}
									elsif ($tmp[$i]=~/6merpair/){delete $tmp[$i];}
								}
								$core_aln_param=join ("\\\-\\\-",@tmp);
							}
							#print "CORE PARAMS:$core_aln_param\n";<STDIN>; # QA
							my $mafftFormatCoreTreeFile=$dir."PRUNE_BP_FOR_CORE_ALN/tree_$countTrees/".$dataset.".".$orig_prog.".semphy.tree_".$countTrees."CORE.mafftFormat";
							my $mafftCoreMSA=$dir."tree_$countTrees/".$dataset.".".$orig_prog.".tree_$countTrees.".$prog.".aln.CORE";
							$cmd="$aln_program $core_aln_param --retree 1 --treein $mafftFormatCoreTreeFile $aminoFile > $mafftCoreMSA; ";
							# use the core MSA and full tree to build the full alignment
							$cmd.="$aln_program $align_param --retree 1 --treein $mafftFormatTreeFile $mafftCoreMSA > $alnFile";
						}
						else
						{
							$cmd = "($aln_program $align_param --retree 1 --treein $mafftFormatTreeFile $aminoFile > $alnFile)";# >& $stdFile";
						}
					} elsif (uc($prog) eq "CLUSTALW") {
						$cmd = "$aln_program $align_param -infile=$aminoFile -usetree=$treeFile -outfile=$alnFile > $stdFile";
					} elsif (uc($prog) eq "MUSCLE") {
						my $rooted_tree=$treeFile.".rooted";
						unless (-e $rooted_tree) {
							return ("ERROR: file does not exist: $rooted_tree\n");
						}
						$cmd = "$aln_program $align_param -in $aminoFile -usetree_nowarn $rooted_tree -out $alnFile > $stdFile";
					} elsif (uc ($prog) eq "PAGAN") {
						my $rooted_tree=$treeFile.".rooted";
						unless (-e $rooted_tree) {
							return ("ERROR: file does not exist: $rooted_tree\n");
						}
						$cmd = "$aln_program $align_param --seqfile $aminoFile --treefile $rooted_tree --outfile $alnFile > $stdFile";
					}
					#print "$cmd\n"; # QA
					system ("$cmd");
					if ($update_file ne "") {
					    open (PROGRESS,">$update_file");
					    print PROGRESS "\n<ul><li>",$countTrees+1," out of $bp_repeats alternative alignments were created</li></ul>\n";
					    close (PROGRESS);
					}
#			print "$cmdFile\n";
#			open IN, ">$cmdFile" or return ("cannot open file $cmdFile\n");
#			print IN "#!/bin/sh\n\ncd ".$dir."tree_$countTrees\n\n";
#			print IN "$cmd\n";
#			close IN;
#			system ($cmd);
					if ($prog eq "PRANK") {
						if ($PRANK_VERSION==0)
						{
							system ("cp $alnFile.1.fas $alnFile");
						}
						elsif ($PRANK_VERSION eq "121218")
						{
							system ("cp $alnFile.best.fas $alnFile");
						} 
						else # ver >=120626
						{
							system ("cp $alnFile.2.fas $alnFile");
						}
					}
					if ($prog eq "CLUSTALW")
					{
						MSA_parser::convert_msa_format($alnFile,"clustalw","$alnFile.fs","fasta");
						system ("mv $alnFile $alnFile.orig");
						system ("mv $alnFile.fs $alnFile");
					}
					if ($prog eq "PAGAN")
					{
						move("$alnFile".".fas", "$alnFile");
					}
					#my $alias = $dataset.".".$prog;
					#my $qsub = "qsub -q heavy -N $alias $cmdFile";
					# print $qsub."\n";
					#system ($qsub);
				}
			}
			exit 0;
		} else {
			die "ERROR: fork failed: $!\n";
		}

	}

	# Wait for child processes to end
	foreach (@children) {
        my $pid = waitpid($_, 0);
		print "done with pid $pid\n";
		#if ($update_file ne "") {
 		#	open (PROGRESS,">$update_file");
 		#	print PROGRESS "\n<ul><li>",$countTrees+1," out of $bp_repeats alternative alignments were created</li></ul>\n";
 		#	close (PROGRESS);
		#}
		
	}

	return "ok";
}


sub rootTree {
####################################################################################################################
#die "USAGE: $0 inTree outTree
#inTree must be an unrooted tree - i.e. root node has at least 3 sons.
#In outTree the root will have 2 sons
#(all direct sons of root will be made biforcating - the rest of tree is left untouched)\n"
####################################################################################################################
	my $inTree=shift;
	my $outTree=shift;
	
	my $in  = new Bio::TreeIO(-file   => "$inTree",
                                  -format => "newick");
	my $out = new Bio::TreeIO(-file   => ">$outTree",
                                  -format => "newick");

	while (my $tree = $in->next_tree) {
		my $root = $tree->get_root_node;
		my @sons = $root->each_Descendent;

		# Remove edges between root-sons
		foreach my $son (@sons) {
			$root->remove_Descendent($son);
		}

        # Iteratively add 
		my $currFather = $root;
		while (@sons > 2) {
			my $son = shift @sons;
			$currFather->add_Descendent($son);
			my $midNode = new Bio::Tree::Node();
			$currFather->add_Descendent($midNode);
			$midNode->branch_length(0);
			$currFather = $midNode;
		}
		$currFather->add_Descendent($sons[0]);
		$currFather->add_Descendent($sons[1]);
		$sons[0]->branch_length(0);
		$sons[1]->branch_length(0);

		$out->write_tree($tree);
	}
}

sub pullOutBPtrees {
####################################################################################################################
# pull out all the BP NJ trees into the BP directory
# pull out the original NJ tree (that was done on the complete MSA file)
####################################################################################################################

#die "Usage: $0 <noBPdir> <dataset> <num_BP_repeats> <alignment program>" if (scalar(@ARGV) < 4);

	my $noBPdir = shift;
	my $dataset = shift; 
	my $bp_repeats = shift;
	my $alnProg = shift;

	unless ($noBPdir =~ m/\/$/) {
    	$noBPdir.="/";
	}

	my $BPdir = $noBPdir."BP/";
	unless (-e  $BPdir) {
    	system ("mkdir $BPdir");
	}
	my $nonUniqueTreesDir="";
	if ($alnProg ne "MAFFT") # BUILED ALIGNEMT ONLY FOR UNIQUE TREES
	{	
		$nonUniqueTreesDir = $BPdir."nonUniqueTrees/";
		unless (-e  $nonUniqueTreesDir) {
			system ("mkdir $nonUniqueTreesDir");
		}
	}

#	my $robinsonFouldCommand = "$phylonet_prog rf";
	# open semphy file
	my $semphyLogFile = $BPdir.$dataset.".".$alnProg.".semphy.log";
	print "semphy log file: $semphyLogFile\n";
	open IN, "<$semphyLogFile" or return "can't open file $semphyLogFile";
	my $countTrees=0;
	my $countUniqueTrees=0;
	my @numRepeats;
	my $readTree=0; # this is a flag. 2 == read the tree. 3 == read the tree and check for uniqueness.
	my $treeDir="";
	my $treeFile="";
	foreach my $line (<IN>) {
		if ($line =~ m/^\# Finished tree reconstruction\./) {
			$readTree=1;
		} 
		elsif (($readTree==1) && ($line =~ m/^\# The tree/)) {
			$treeFile = $noBPdir.$dataset.".".$alnProg.".semphy.tree";
			$readTree=2;
		}
		elsif (($readTree==1) && ($line =~ m/The reconsructed tree/)) {
			$readTree=2;
			if ($alnProg eq "MAFFT")
			{
				$treeDir = $BPdir."/tree_".$countTrees."/"
			}
			else # CHECK UNIQUE TREE ONLY NOT FOR MAFFT
			{
				$treeDir = $nonUniqueTreesDir."/tree_".$countTrees."/";
			}
			unless (-e  $treeDir) {system ("mkdir $treeDir");}
			$treeFile = $treeDir.$dataset.".".$alnProg.".semphy.tree_".$countTrees;
			$countTrees++;
			$readTree=3 if ($alnProg ne "MAFFT"); # CHECK UNIQUE TREE ONLY NOT FOR MAFFT
		}
		elsif ($readTree>=2) {
			# open treefile
			open OUT, ">$treeFile" or return "can't open file $treeFile";
			# write the tree into the treefile
			print OUT $line;
			# close treefile
			close OUT;
			if ($readTree==3) { # compare the tree to all other unique trees
				for ($i=0;$i<$countUniqueTrees;++$i) {
					$uniqueTreeFile = $BPdir."tree_".$i."/".$dataset.".".$alnProg.".semphy.tree_".$i;
#					$RobinsonFouldResFile = $treeDir."rf.".$i.".std";
					$isEqualTopologyResFile = $treeDir."isEqualTopology.".$i.".std";
#					my $rfCommand = "$robinsonFouldCommand -m $treeFile -e $uniqueTreeFile";
					my $isEqualTopologyCommand = "$isEqualTopologyProg $treeFile $uniqueTreeFile";
#					my @rfResults = `$rfCommand`;
					my $isEqualTopology = `$isEqualTopologyCommand`;
#					my $editedRFresults = $rfResults[0];
#					open OUT_RF, ">$RobinsonFouldResFile" or return "can't open file $RobinsonFouldResFile";
#					print OUT_RF "@rfResults";
					open OUT_EQUAL_TOP, ">$isEqualTopologyResFile" or return "can't open file $isEqualTopologyResFile";
					print OUT_EQUAL_TOP "$isEqualTopology\n";
					close OUT_EQUAL_TOP;
#					if ($editedRFresults =~ m/ERROR/) {
#						print "skipping error in rfResults of $treeFile and $uniqueTreeFile : $editedRFresults\n";
#						print OUT_RF "skipping error in rfResults of $treeFile and $uniqueTreeFile : $editedRFresults\n";
#						close OUT_RF;
#						next;
#					}
#					close OUT_RF;
#					print "isEqualTopology == $isEqualTopology\n"; #debug OP
#					if ($editedRFresults =~ m/^0\s+0\s+/) { # same tree, FP=0 and FN=0
					if ($isEqualTopology == 1) { # same tree
						$numRepeats[$i]++;
#						print "$isEqualTopology == 1\n"; # debug OP
						last;
					}
					if ($isEqualTopology == 2) {
						print "skipping ERROR in isEqualTopology of $treeFile and $uniqueTreeFile\n";
						next;
					}
				}
				if ($i == $countUniqueTrees) { # The new tree is unique
					push (@numRepeats,1);
					my $uniqueTreeDir = $BPdir."tree_".$countUniqueTrees."/";
					unless (-e  $uniqueTreeDir) {system ("mkdir $uniqueTreeDir");}
					my $uniqueTreeFile = $uniqueTreeDir.$dataset.".".$alnProg.".semphy.tree_".$countUniqueTrees;
					system ("cp $treeFile $uniqueTreeFile");
					$countUniqueTrees++;
				}
			}
			$readTree=1;
		}
	}
	close IN;
	if ($countTrees != $bp_repeats) {
		return "ERROR: dataset: $dataset \t countTrees: $countTrees while it should be $bp_repeats \n";
	}
	if ($alnProg ne "MAFFT")
	{
		# print the number of repeats per unique tree into file
		my $numRepeatsFile = $BPdir."numRepeats";
		open OUT_NUM_REPEATS, ">$numRepeatsFile" or return "can't open file $numRepeatsFile";
		print OUT_NUM_REPEATS "@numRepeats";
		close OUT_NUM_REPEATS;
		return "ok",$countUniqueTrees,\@numRepeats;
	}
	else
	{
		return "ok";	
	}
}

sub copyBootstrapMSA2oneDir {

#die "Usage: $0 <BP_MSA_dir> <BPdir> <alignment program (prank/mafft/clustal/muscle) <ref2array - numRepeats4UniqueTree>" if (scalar(@ARGV) < 4);

	my ($BP_MSA_dir,$BPdir,$prog,$RefNumRepeats4UniqueTree) = @_;
	my @numRepeats4UniqueTree= @$RefNumRepeats4UniqueTree;
	print "numRepeats4UniqueTree: @numRepeats4UniqueTree\n";
	mkdir ($BP_MSA_dir);
	my $countTrees=0;
	for ($i=0; $i<scalar(@numRepeats4UniqueTree); ++$i) {
		$numRepeats = $numRepeats4UniqueTree[$i];
		for (my $j=0; $j < $numRepeats; ++$j) {
			system ("cp $BPdir/tree_$i/*.$prog.aln $BP_MSA_dir/repeated_tree_$countTrees.$prog.aln");
			$countTrees++;
		}
	}

	return "ok";
}
	
sub codes2nameFastaFrom1 {
	my $Aln_with_Codes=shift;;
	my $codes_file=shift;
	my $Aln_with_names=shift;
	my %Codes=();
	open (CODES,$codes_file) or return ("Guidance::codes2nameFastaFrom1 Can't open the Codes file: '$codes_file' $!");
	while (my $line=<CODES>)
	{
		chomp $line;
		my ($Seq_name,$Code)=split("\t",$line);
		$Codes{$Code}=$Seq_name;
	}
	close (CODES);
	open (IN,$Aln_with_Codes) or return ("Guidance::codes2nameFastaFrom1 Can't open the Alignment with codes file: '$Aln_with_Codes' $!");
	open (OUT,">$Aln_with_names") or return ("Guidance::codes2nameFastaFrom1 Can't open the Alignment with names file: '$Aln_with_names' for writing $!");
	while (my $line=<IN>){
		chomp ($line);
		if ($line=~/^>([0-9]+_placeholder_to_prevent_PRANK_error)/) #RvdL
		{
			print OUT ">",$Codes{$1},"\n";
		}
		else # Seq Line
		{
			print OUT "$line\n";
		}
	}
	close OUT;
	close IN;
	return "OK";
}	
	

sub ConvertNamesOfAlignWithSeed
###################################################################################################
# Take alignment created with seed (by MAFFT) and convert the seq with _seed_ prefix into numbers
# IMPORTANT: The seed seq must be first so MAFFT --reorder argument must not be use
###################################################################################################
{
	my $aln=shift;
	my $out=shift;
	open (IN,$aln) || return  "Guidance::ConvertNamesOfAlignWithSeed: Can't open In '$aln' $!";
	open (OUT,">$out") || return "Guidance::ConvertNamesOfAlignWithSeed: Can't open Out '$out' $!";
	my $seed_counter=0;
	while (my $line=<IN>)
	{
		if ($line=~/>_seed_(.*)/)
		{
			$seed_counter++;
			$line=~s/>_seed_(.*)/>$seed_counter/;
			print OUT "$line";
		}
		else
		{
			print OUT $line;
		}
	}
	close (IN);
	close (OUT);
	return "ok";
}
sub name2codeFasta_without_codded_out {
####################################################################################################################
# Convert the names in a fasta file to numbers, and creates a code file with the names and the codes (running number)
###################################################################################################################
	my $in_fileName = shift;
	my $code_fileName = shift;
	my $counter_offset= shift; # optional

	my $in_file = Bio::SeqIO->new(-file => $in_fileName , '-format' => 'Fasta');
	my $code_file = new FileHandle(">>$code_fileName") or return ("Can't write to $code_fileName $!");
	$counter_offset=1 if (!defined $counter_offset);
	$counter_offset=1 if ($counter_offset==0);
	my $counter = $counter_offset;
	my $i;

	while ( my $seqObj = $in_file->next_seq() ) {
		my $name = $seqObj->display_id();
		$name.= " ".$seqObj->desc()   if ($seqObj->desc());
		print $code_file "$name\t$counter\n";
		$counter++;
	}
	$in_file->close();
	$code_file->close();
	return "ok",$counter;
}	
sub name2codeFastaFrom1 {
####################################################################################################################
# Convert the names in a fasta file to numbers, and creates a code file with the names and the codes (running number)
###################################################################################################################
	my $in_fileName = shift;
	my $code_fileName = shift;
	my $out_fileName = shift;
	my $counter_offset=shift; # optional

	my $in_file = Bio::SeqIO->new(-file => $in_fileName , '-format' => 'Fasta');
	my $code_file = new FileHandle(">>$code_fileName") or return ("Can't write to $code_fileName $!");
	my $out_file = new FileHandle(">$out_fileName") or return ("Can't write to $out_fileName");
	$counter_offset=1 if (!defined $counter_offset);
	$counter_offset=1 if ($counter_offset==0);
	my $counter = $counter_offset;
	my $i;

	while ( my $seqObj = $in_file->next_seq() ) {
		my $name = $seqObj->display_id();
    		$name.= " ".$seqObj->desc()   if ($seqObj->desc());
    		print $code_file "$name\t$counter" . "_placeholder_to_prevent_PRANK_error\n"; #RvdL
    		my $seq = $seqObj->seq();
    		print $out_file ">$counter" . "_placeholder_to_prevent_PRANK_error\n"; #RvdL
        	for($i=0;$i<length($seq);$i+=60){
                	print $out_file substr($seq,$i,60) . "\n";
        	}
        	if($i<length($seq)){
                	print $out_file substr($seq,$i,length($seq)-$i);
        	}
        	print $out_file "\n";
    		$counter++;
	}
	$out_file->close();
	$in_file->close();
	$code_file->close();
	return "ok";
}

sub average{
        my $data  = $_[0];
		my $data_size=scalar(@{$data});
        if ($data_size==0) {
			return ("ERROR: Empty array\n");
        }
        my $total = 0;
        foreach my $val (@{$data}) {
			$total += $val;
        }
        my $average = $total / $data_size;
        return $average;
}
sub stdev{
        my $data  = $_[0];
		my $average = $_[1];
		my $data_size=scalar (@{$data});
        if($data_size == 1){
			return 0;
        }
        my $sqtotal = 0;
        foreach my $val (@{$data}) {
			$sqtotal += ($average-$val) ** 2;
        }
        my $std = ($sqtotal / ($data_size-1)) ** 0.5;
        return $std;
}
sub Calculate_mean_and_std
# get a delimited file and calcilate the std and average
{
	my $DataFile=shift;
	my $Column=shift;
	$Column--; # for split 
	my $AVERAGE="NAN";
	my $STD="NAN";
	my @data=();
	open (DATA_FILE,"<$DataFile") or return ("Calculate_mean_and_std: Can't open file '$DataFile' $!");
	my $line=<DATA_FILE>; # for header
	while ($line=<DATA_FILE>)
	{
		chomp ($line);
		$line=~ s/^\s+|\s+$//g;
		my @line=split (/\s+/,$line);
		push (@data,$line[$Column]) if ((defined $line[$Column]) and ($line[$Column]=~/\d/));
	}
	close (DATA_FILE);
	my $avg=average(\@data);
	my $std=stdev(\@data,$avg);
	return ($avg,$std);
}
sub removeLowSPseq_Consider_Z{
	# will remove all sequences in which their Z score is bellow cutoff and their SP score is also below cutoff
	my $msaFile=shift;
  	my $SeqSpFile=shift;
  	my $outFile=shift;
  	my $cutoof=shift;
	my $Z_cutoff=shift;
 	my $removed_seq_file=shift;

	my ($mean,$std)=Calculate_mean_and_std ($SeqSpFile,2);
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $msaFile) or die "Can't open $msaFile: $!";
	my $aln = $in->next_aln;
	$aln->verbose(1); #HAIM COMMNET
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat(); #HAIM COMMENT
	
  	open (SEQ_SP_SCORES,"<$SeqSpFile") or return ("removeLowSPseq_Consider_Z: Can't open file '$SeqSpFile' $!");
  	open (OUT,">$outFile") or return ("removeLowSPseq_Consider_Z: Can't open file '$outFile' $!");
  	open (REMOVED_SEQ,">$removed_seq_file") or return ("removeLowSPseq_Consider_Z: Can't open file '$removed_seq_file' $!");
	my $line=<SEQ_SP_SCORES>; #For header
	while ($line=<SEQ_SP_SCORES>){
		if ($line=~m/^\s*(\d+)\s+(\d+(\.\d+)?)/)
		{
			my $row_num=$1;
			my $Seq_SP_Score=$2;
			my $seq_obj = $aln->get_seq_by_pos($row_num);
			my $seq=$seq_obj->seq();
			$seq=~s/-//g;
			$seq =~ s/(.{60,60})/$1\n/g ;
  			$seq .= "\n" unless (substr($seq, -1, 1) eq "\n") ;
			my $Z_score="NaN";
			if ($std>0) {$Z_score=($Seq_SP_Score-$mean)/$std;}
			if ($Seq_SP_Score eq "nan")
			{
				print OUT ">".$seq_obj->id(),"_SP_$Seq_SP_Score_Z_$Z_score\n",$seq,"\n";
			}
			elsif (($Z_score ne "NaN") and ($Z_score<$Z_cutoff)) # a negative outlier
			{
				if (($Seq_SP_Score<$cutoof)) # To avoid filter highly scored position
				{
					print REMOVED_SEQ ">".$seq_obj->id(),"_SP_$Seq_SP_Score"."_Z_$Z_score\n",$seq,"\n";	
				}
				else
				{
					print OUT ">".$seq_obj->id(),"_SP_$Seq_SP_Score"."_Z_$Z_score\n",$seq,"\n";
				}
			}
			else
			{
				print OUT ">".$seq_obj->id(),"_SP_$Seq_SP_Score"."_Z_$Z_score\n",$seq,"\n";
			}
		}
 	}
	close (SEQ_SP_SCORES);
	close (OUT);
	close (REMOVED_SEQ);
	return "OK";

}

sub removeLowSPseq{
  	my $msaFile=shift;
  	my $SeqSpFile=shift;
  	my $outFile=shift;
  	my $cutoof=shift;
 	my $removed_seq_file=shift;
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $msaFile) or die "Can't open $msaFile: $!";
	my $aln = $in->next_aln;
	$aln->verbose(1); #HAIM COMMNET
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat(); #HAIM COMMENT
	
  	open (SEQ_SP_SCORES,"<$SeqSpFile") or return ("removeLowSPseq: Can't open file '$SeqSpFile' $!");
  	open (OUT,">$outFile") or return ("removeLowSPseq: Can't open file '$outFile' $!");
  	open (REMOVED_SEQ,">$removed_seq_file") or return ("removeLowSPseq: Can't open file '$removed_seq_file' $!");
	my $line=<SEQ_SP_SCORES>; #For header
	while ($line=<SEQ_SP_SCORES>){
		if ($line=~m/^\s*(\d+)\s+(\d+(\.\d+)?)/)
		{
			my $row_num=$1;
			my $Seq_SP_Score=$2;
			my $seq_obj = $aln->get_seq_by_pos($row_num);
			my $seq=$seq_obj->seq();
			$seq=~s/-//g;
			$seq =~ s/(.{60,60})/$1\n/g ;
  			$seq .= "\n" unless (substr($seq, -1, 1) eq "\n") ;
			if ($Seq_SP_Score>=$cutoof)
			{
				print OUT ">".$seq_obj->id(),"\n",$seq,"\n";
			}
			else
			{
				print REMOVED_SEQ ">".$seq_obj->id(),"\n",$seq,"\n";	
			}
		}
 	}
	close (SEQ_SP_SCORES);
	close (OUT);
	close (REMOVED_SEQ);
	return "OK";
  }

sub removeLowSPsites_Consider_Z {
	# no checks of the input files are done.
	#die "USAGE: $0 MSA_FILE SP_FILE OUT_FILE CUTOFF Z_CUTOFF" if (@ARGV < 5);

	my $msaFile=shift;
	my $spFile=shift;
	my $outFile=shift;
	my $cutoff=shift;
	my $Z_cutoff=shift;
	my $Pos_removed_file=shift;
	  
	my $MSA_Length=0;
	my $Num_Pos_removed=0;
	my $in_fasta = Bio::SeqIO->new(-file => $msaFile, '-format' => 'fasta');
	my @seqs;

	# read the file into an array
	while (my $seqObj = $in_fasta->next_seq()) {
    		push(@seqs,$seqObj);
		if ($MSA_Length==0)
		  {
		    $MSA_Length=$seqObj->length;
		  }
	}
	my ($mean,$std)=Calculate_mean_and_std ($spFile,2);
	open IN, "<$spFile" or return ("removeLowSPsites_Consider_Z: can't open file '$spFile'\n");
	open OUT, ">$outFile" or return ("removeLowSPsites_Consider_Z: can't open file '$outFile' $!\n");
	if ($Pos_removed_file ne ""){open REMOVED_POS, ">$Pos_removed_file" or return ("removeLowSPsites_Consider_Z: can't open file '$Pos_removed_file' $!\n");}
	my $numRemovedPos = 0;
	foreach my $line(<IN>) {
		if ($line =~ m/^\s*(\d+)\s+(\d+(\.\d+)?)/) {
			my $site_num=$1;
			my $site_score=$2;
			my $Z_score="NaN";
			if ($std>0) {$Z_score=($site_score-$mean)/$std;}
			if (($Z_score ne "NaN") and ($Z_score<-$Z_cutoff)) # an outlier
			{
				if ($site_score < $cutoff) {
					removePos($site_num,$numRemovedPos,\@seqs);
					print REMOVED_POS "Remove Pos: $site_num\tScore: $site_score\tZ_Score:$Z_score\n" if ($Pos_removed_file ne "");
					$numRemovedPos++;
				}
			}
		}
	}
	foreach my $seqObj (@seqs) { 
    		my $id = $seqObj->id();
    		my $seq = $seqObj->seq();
    		print OUT ">$id\n$seq\n";
	}
	close (IN);
	close (OUT);
	close (REMOVED_POS);
	return ("OK",$numRemovedPos,$MSA_Length);
}
sub removeLowSPsites {
	# no checks of the input files are done.
	#die "USAGE: $0 MSA_FILE SP_FILE OUT_FILE CUTOFF" if (@ARGV < 4);

	my $msaFile=shift;
	my $spFile=shift;
	my $outFile=shift;
	my $cutoff=shift;
	my $Pos_removed_file=shift;
#	my $Alphabet=shift;
	my $MSA_Length=0;
	my $Num_Pos_removed=0;
	my $in_fasta = Bio::SeqIO->new(-file => $msaFile, '-format' => 'fasta');
#	my $in_fasta = Bio::SeqIO->new(-file => $msaFile, '-format' => 'fasta',-alphabet => $Alphabet);

	my @seqs;

	# read the file into an array
	while (my $seqObj = $in_fasta->next_seq()) {
		push(@seqs,$seqObj);
		if ($MSA_Length==0)
		{
		    $MSA_Length=$seqObj->length;
		}
	}
	open IN, "<$spFile" or return ("removeLowSPsites: can't open file '$spFile'\n");
	open OUT, ">$outFile" or return ("removeLowSPsites: can't open file '$outFile' $!\n");
	if ($Pos_removed_file ne ""){open REMOVED_POS, ">$Pos_removed_file" or return ("removeLowSPsites: can't open file '$Pos_removed_file' $!\n");}
	
	my $numRemovedPos = 0;
	foreach my $line(<IN>) {
		if ($line =~ m/^\s*(\d+)\s+(\d+(\.\d+)?)/) {
			if ($2 < $cutoff) {
				removePos($1,$numRemovedPos,\@seqs);
#				removePos($1,$numRemovedPos,\@seqs,$Alphabet);
				print REMOVED_POS "Remove Pos: $1\tScore: $2\n" if ($Pos_removed_file ne "");
				$numRemovedPos++;
			}
		}
	}
	foreach my $seqObj (@seqs) { 
		my $id = $seqObj->id();
		my $seq = $seqObj->seq();
		print OUT ">$id\n$seq\n";
	}
	close (IN);
	close (OUT);
	close (REMOVED_POS);
	return ("OK",$numRemovedPos,$MSA_Length);
}
sub removePos {
    my $numRemovedPos=$_[1];
    my $seq_ref=$_[2];
#	my $Alphabet=$_[3];
    my $pos2remove = $_[0] - $numRemovedPos;
	foreach my $seqObj (@{$seq_ref}){
		my $new_seq = "";
		if ($pos2remove>1) {
			$new_seq = $seqObj->subseq(1,$pos2remove-1);
		}
		if ($pos2remove< $seqObj->length()) {
			$new_seq .= $seqObj->subseq($pos2remove+1,$seqObj->length());
		}
		$seqObj->seq($new_seq);
#		$seqObj->alphabet($Alphabet);
    } 
}

sub removeLowSPsites_NoBioPerl_Consider_Z {
	# no checks of the input files are done.
	#die "USAGE: $0 MSA_FILE SP_FILE OUT_FILE CUTOFF Z_CUTOFF POS_REMOVED_FILE" if (@ARGV < 6);

	my $msaFile=shift;
	my $spFile=shift;
	my $outFile=shift;
	my $cutoff=shift;
	my $Z_cutoff=shift;
	my $Pos_removed_file=shift;
	  
	my $Num_Pos_removed=0;
	my @ans=readMSA_to_Hash($msaFile);
	if ($ans[0] ne "OK"){return "removeLowSPsites_Consider_Z:$ans[0]\n";}
	my $MSA_HashRef=$ans[1]; # hash of array, each key is seq ID and the value is an array with the seq
	my $MSA_Length=$ans[2];
	my $MSA_Order_ArrayRef=$ans[3];
	my ($mean,$std)=Calculate_mean_and_std ($spFile,2);
	open IN, "<$spFile" or return ("removeLowSPsites_Consider_Z: can't open file '$spFile'\n");
	open OUT, ">$outFile" or return ("removeLowSPsites_Consider_Z: can't open file '$outFile' $!\n");
	if ($Pos_removed_file ne ""){open REMOVED_POS, ">$Pos_removed_file" or return ("removeLowSPsites_Consider_Z: can't open file '$Pos_removed_file' $!\n");}
	my $numRemovedPos = 0;
	foreach my $line(<IN>) {
		if ($line =~ m/^\s*(\d+)\s+(\d+(\.\d+)?)/) {
			my $site_num=$1;
			my $site_score=$2;
			my $Z_score="NaN";
			if ($std>0) {$Z_score=($site_score-$mean)/$std;}
			if (($Z_score ne "NaN") and ($Z_score<$Z_cutoff)) # an outlier
			{
				if ($site_score < $cutoff) {
					$MSA_HashRef=removePos_noBioPerl($site_num-1,$MSA_HashRef);
					print REMOVED_POS "Remove Pos: $site_num\tScore: $site_score\tZ_Score: $Z_score\n" if ($Pos_removed_file ne "");
					$numRemovedPos++;
				}
			}
		}
	}
	foreach my $seqID (@{$MSA_Order_ArrayRef}) { 
		my $seq = join ("",@{$MSA_HashRef->{$seqID}});
		if ($seq=~/^[-]+$/) {print "WARNNING: After removing positions scored below $cutoff and Z_Score: $Z_Score, the sequence $seqID comprised of only gap characters...\n";}
		print OUT ">$seqID\n$seq\n";
	}
	close (IN);
	close (OUT);
	close (REMOVED_POS);
	return ("OK",$numRemovedPos,$MSA_Length);
}
sub removeLowSPsites_NoBioPerl {
	# no checks of the input files are done.
	#die "USAGE: $0 MSA_FILE SP_FILE OUT_FILE CUTOFF POS_REMOVED_FILE" if (@ARGV < 5);

	my $msaFile=shift;
	my $spFile=shift;
	my $outFile=shift;
	my $cutoff=shift;
	my $Pos_removed_file=shift;
	
	my $Num_Pos_removed=0;
	my @ans=readMSA_to_Hash($msaFile);
	if ($ans[0] ne "OK"){return "removeLowSPsites:$ans[0]\n";}
	my $MSA_HashRef=$ans[1];
	my $MSA_Length=$ans[2];
	my $MSA_Order_ArrayRef=$ans[3];
	open IN, "<$spFile" or return ("removeLowSPsites: can't open file '$spFile'\n");
	open OUT, ">$outFile" or return ("removeLowSPsites: can't open file '$outFile' $!\n");
	if ($Pos_removed_file ne ""){open REMOVED_POS, ">$Pos_removed_file" or return ("removeLowSPsites: can't open file '$Pos_removed_file' $!\n");}
	
	my $numRemovedPos = 0;
	foreach my $line(<IN>) {
		chomp ($line);
		if ($line =~ m/^\s*(\d+)\s+(\d+(\.\d+)?)/) {
			if ($2 < $cutoff) {
				my $Pos=$1;
				$MSA_HashRef=removePos_noBioPerl ($Pos-1,$MSA_HashRef); # Pos-1 because array start from 0;
				print REMOVED_POS "Remove Pos: $1\tScore: $2\n" if ($Pos_removed_file ne "");
				$numRemovedPos++;
			}
		}
	}
	foreach my $seqID (@{$MSA_Order_ArrayRef}) { 
		my $seq = join ("",@{$MSA_HashRef->{$seqID}});
		if ($seq=~/^[-]+$/) {print "WARNNING: After removing positions scored below $cutoff, the sequence $seqID comprised of only gap characters...\n";}
		print OUT ">$seqID\n$seq\n";
	}
	close (IN);
	close (OUT);
	close (REMOVED_POS);
	return ("OK",$numRemovedPos,$MSA_Length);
}
sub printMSA_Hash
{
	my $refToHash=shift;
	foreach my $key (keys %{$refToHash})
	{
		print ">$key\n",join ("",@{$refToHash->{$key}}),"\n";
	}
}
sub readMSA_to_Hash
{
	# Take an MSA in FASTA format and return (1) an hash of arrays where seq id is key and the seq is an array; (2) MSA length (3) ref to array with the order of seq header in orig file
	my $MSA_File=$_[0];
	my %MSA_Hash=();
	my $MSA_Length=0;
	my @MSA_Order=();
	# open file
	open (my $inMSA, "<", $MSA_File)  or return ("GUIDANCE::readMSA_to_Hash: FAILED to open '$MSA_File' $!");
	## 1.1. Read FASTA header and save it
	my $fastaLine = <$inMSA>;
	while (defined $fastaLine) {
		chomp $fastaLine;
		my $header = substr($fastaLine,1);
		## 1.2. Read seq until next header
		$fastaLine = <$inMSA>;
		my $seq = "";
		while ((defined $fastaLine) and
			   (substr($fastaLine,0,1) ne ">" )) {
			chomp $fastaLine;
			$seq .= $fastaLine;
			$fastaLine = <$inMSA>;
		}
		## 2.1 update hash
		my @seq_arr=split(//,$seq);
		$MSA_Hash{$header}=[@seq_arr];
		push (@MSA_Order,$header);
		if ($MSA_Length==0)
		{
			$MSA_Length=length ($seq);
		}
	}
	# close file
	close ($inMSA); 
	return ("OK",\%MSA_Hash,$MSA_Length,\@MSA_Order);
}

sub removePos_noBioPerl {
# Get MSAhash where each key is seqID and each value is reff to array with seq;
# Give a position, rumove the specific char from all sequences
	my $PosToRemove=$_[0];
	my $MSA_Hash_ref=$_[1];
	foreach my $Seq_ID (keys %{$MSA_Hash_ref})
	{
		$MSA_Hash_ref->{$Seq_ID}[$PosToRemove]="";
	}
	return $MSA_Hash_ref;
}

sub codes2name_scoresFile
{
	my $Score_File=shift;
	my $Codes_File=shift;
	my $Out=shift;
	
	my %Codes=();
	open (CODES,$Codes_File) or return ("codes2name_scoresFile:Can't open the Codes file: '$Codes_File' $!");
	while (my $line=<CODES>)
	{
		chomp $line;
		my ($Seq_name,$Code)=split("\t",$line);
		$Codes{$Code}=$Seq_name;
	}
	close (CODES);
	open (OUT,">$Out") or return ("codes2name_scoresFile: Can't open out file: '$Out' $!");
	open (SCORES,$Score_File) or return ("codes2name_scoresFile: Can't open Scores file: '$Score_File' $!");
	my $line=<SCORES>; #Header;
	print OUT "SEQUENCE_NAME\tSEQUENCE_SCORE\n";
	while ($line=<SCORES>)
	{
		my $code="";
		my $score="";
		if ($line=~/([0-9]+)\s+([0-9.]+)/){
			$code=$1;
			$score=$2;
		}
		if ($Codes{$code})
		{
			print OUT "$Codes{$code}\t$score\n";
		}
		else
		{
			print OUT "$line";
		}
	}
	close (OUT);
	close (SCORES);
	
	return "OK";
}
  
sub Convert_to_Codons_Numbering
{
	my $score_file=shift;
	my $score_codons_file=shift;
	
	open (IN,$score_file)|| return ("Can't open the input file: '$score_file' $!\n");
	open (OUT,">$score_codons_file") || return ("Can't open the output file: '$score_codons_file' $!\n");
	while (my $line=<IN>)
	{
		chomp($line);
		$line=trim($line);
		my @line=split(/\s+/,$line);
#	  my $array_size=scalar(@line);
		if ($line[0]=~/[0-9]+/)
		{
			$array_size=scalar(@line)-1;
			print OUT $line[0]*3-2,"\t",join("\t",@line[1..$array_size]),"\n";
			print OUT $line[0]*3-1,"\t",join("\t",@line[1..$array_size]),"\n";
			print OUT $line[0]*3,"\t",join("\t",@line[1..$array_size]),"\n";
		}
		else 
		{
			print OUT join("\t",@line),"\n";
		}
	}
	close (IN);
	close (OUT);
	return "OK";
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^[\s\t]+//;
	$string =~ s/[\s\t]+$//;
	return $string;
}
sub MSA_row_Num_to_Seq_Name
{
	my $MSA=shift;
	my $Codes_Name_File=shift;
	my %MSA_Row_Seq_Name=();
	my %Code_Names=();
	open (CODES,$Codes_Name_File) or return ("Guidance::MSA_row_Num_to_Seq_Name Can't open the Codes file: '$Codes_Name_File' $!");
	while (my $line=<CODES>)
	{
		chomp $line;
		my ($Seq_name,$Code)=split("\t",$line);
		$Code_Names{$Code}=$Seq_name;
	}
	close (CODES);
	# Read MSA
	my $MSA_Depth;
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $MSA) or return ("MSA_row_Num_to_Seq_Name can't open $inMsaFile: $!");
	my $aln = $in->next_aln;
	$aln->verbose(1);
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat();
	@ans=MSA_parser::check_msa_licit_and_size($MSA,"fasta","no");
	if ($ans[0] eq "OK"){$MSA_Depth=$ans[1];}
	else {return "MSA_row_Num_to_Seq_Name: ".join (" ",@ans);}
	for(my $i=1;$i<=$MSA_Depth;$i++){
		my $seq = $aln->get_seq_by_pos($i);
		my $Seq_Name=$Code_Names{$seq->id};
		$MSA_row_Num_to_Seq_Name{$i}=$Seq_Name;
		print "ROW $i\t$Seq_Name\n";
	}
	return ("OK",\%MSA_row_Num_to_Seq_Name);
}
sub codes2name_scoresFile_NEW
# The Scores file
{
	my $Score_File=shift;
	my $Codes_File=shift;
	my $MSA_File=shift;
	my $Out=shift;

	my %MSA_row_Num_to_Seq_Name=();
	my %Code_Names=();
	
	# Read codes 
	open (CODES,$Codes_File) or return ("Guidance::codes2name_scoresFile_NEW Can't open the Codes file: '$Codes_File' $!");
	while (my $line=<CODES>)
	{
		chomp $line;
		my ($Seq_name,$Code)=split("\t",$line);
		$Code_Names{$Code}=$Seq_name;
	}
	close (CODES);
	# Read MSA to see which seq in which row and assign the correct code name
	my $MSA_Depth;
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $MSA_File) or return ("Guidance::codes2name_scoresFile_NEW can't open $MSA_File: $!");
	my $aln = $in->next_aln;
	$aln->verbose(1);
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat();
	@ans=MSA_parser::check_msa_licit_and_size($MSA_File,"fasta","no");
	if ($ans[0] eq "OK"){$MSA_Depth=$ans[1];}
	else {return ": ".join (" ",@ans);}
	for(my $i=1;$i<=$MSA_Depth;$i++){
		my $seq = $aln->get_seq_by_pos($i);
		my $Seq_Name=$Code_Names{$seq->id};
		$MSA_row_Num_to_Seq_Name{$i}=$Seq_Name;
#		print "ROW $i\t$Seq_Name\n";
	}
	#Add names to score file
	open (OUT,">$Out") or return ("Guidance::codes2name_scoresFile_NEW: Can't open out file: '$Out' $!");
	open (SCORES,$Score_File) or return ("Guidance::codes2name_scoresFile_NEW: Can't open Scores file: '$Score_File' $!");
	my $line=<SCORES>; #Header;
	print OUT "SEQUENCE_NAME\tSEQUENCE_SCORE\n";
	while ($line=<SCORES>)
	{
		my $code="";
		my $score="";
		my $MSA_row_Num="";
		if ($line=~/([0-9]+)\s+([0-9.]+)/){
			$MSA_row_Num=$1;
			$score=$2;
		}
		if (defined $MSA_row_Num_to_Seq_Name{$MSA_row_Num})
		{
			print OUT "$MSA_row_Num_to_Seq_Name{$MSA_row_Num}\t$score\n";
		}
		else
		{
			print OUT "$line";
		}
	}
	close (OUT);
	close (SCORES);
	return "OK";
}

# JALVIEW OUTPUTS
sub make_Jalview_AnnotationGraph
{
	my $Jalview_AnnotFile=shift; # OUT FILE
	my $Data_File=shift;         # file ordered by the X's, first field must be X
	my $Y_label=shift;           # The Y label
	my $Y_data_Col=1;            # The Col of Y data
	my $last_x = 0;
	open (OUT,">$Jalview_AnnotFile") || return ("Can't open outAnnotationsFile '$Jalview_AnnotFile': $!");
	print OUT "JALVIEW_ANNOTATION\n";
	print OUT "BAR_GRAPH\t$Y_label\t";
	
	open (DATA_FILE,$Data_File) || die ("make_Jalview_AnnotationGraph: Can't open data file '$Data_File' $!");
	my $line=<DATA_FILE>; # header
	while ($line=<DATA_FILE>)
	{
		chomp ($line);
		my @data=split(",",$line);
		if (($data[$Y_data_Col] ne "nan") and (($data[0]-$last_x)==1))
		{
			print OUT "$data[$Y_data_Col],$data[$Y_data_Col],$data[$Y_data_Col]|";
			$last_x=$data[0];
		}
		elsif ($data[1] ne "nan") 
		{
			while ($data[0]-$last_x!=1)
			{
				print OUT "0,0,0|";
				$last_x++;
			}
			print OUT "$data[$Y_data_Col],$data[$Y_data_Col],$data[$Y_data_Col]|";
			$last_x=$data[0];
		}
		elsif ($data[$Y_data_Col] eq "nan")
		{
			print OUT "$data[$Y_data_Col],$data[$Y_data_Col],$data[$Y_data_Col]|";
			$last_x=$data[0];
		}
	}
	print OUT "\n";
	close (OUT);
}
sub make_Jalview_Color_MSA
{
	# Data for MSA coloring
	my $inMsaFile=shift;                     # MSA File
	my $scoresFile=shift;                # Scores File
	my $outJalviewFeaturesFile=shift;
	my $codesFile=shift;                 # OPTIONAL
	# Global VARS	
	
	my $sequenceLengthForDisplay=400000;
	
	# Print HTML start
	
	open JALVIEW_FEATURES, ">$outJalviewFeaturesFile" or return ("Can't open outFeaturesFile '$outJalviewFeaturesFile': $!");

#td.Score9{        color: #FFFFFF;        background: #A02560;
#td.Score8{        background: #F07DAB;
#td.Score7{        background: #FAC9DE;
#td.Score6{        background: #FCEDF4;
#td.Score5{        background: #FFFFFF;
#td.Score4{        background: #EAFFFF;
#td.Score3{        background: #D7FFFF;
#td.Score2{        background: #8CFFFF;
#td.Score1{        background: #10C8D1;
#td.ScoreNaN	    background: #C0C0C0;
	print JALVIEW_FEATURES "Score1\t10C8D1\n";
	print JALVIEW_FEATURES "Score2\t8CFFFF\n";
	print JALVIEW_FEATURES "Score3\tD7FFFF\n";
	print JALVIEW_FEATURES "Score4\tEAFFFF\n";
	print JALVIEW_FEATURES "Score5\tFFFFFF\n";
	print JALVIEW_FEATURES "Score6\tFCEDF4\n";
	print JALVIEW_FEATURES "Score7\tFAC9DE\n";
	print JALVIEW_FEATURES "Score8\tF07DAB\n";
	print JALVIEW_FEATURES "Score9\tA02560\n";
	print JALVIEW_FEATURES "ScoreNaN\tC0C0C0\n";
	
	print JALVIEW_FEATURES "STARTGROUP\tGUIDANCE\n";
	
	open SCORES, $scoresFile or return ("Can't open $scoresFile: $!");
	my %scores;
	my %Code_Names;
	foreach (<SCORES>) {
		next if (/^\#/);
		s/^\s+//;
		my ($col, $seq, $score) = split;
		$scores{$seq}[$col] = $score;
	}
	# Read Codes
	if ($codesFile ne "")
	{
		open (CODES,$codesFile) or return ("Guidance::printColoredAlignment Can't open the Codes file: '$codesFile' $!");
		while (my $line=<CODES>)
		{
			chomp $line;
			my ($Seq_name,$Code)=split("\t",$line);
			$Code_Names{$Code}=$Seq_name;
		}
		close (CODES);
	}
	# Read MSA
	my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $inMsaFile) or die "Can't open $inMsaFile: $!";
	my $aln = $in->next_aln;
	$aln->verbose(1); #HAIM COMMNET
	# Otherwise, bioperl adds sequence start/stop values
	$aln->set_displayname_flat(); #HAIM COMMENT
	@ans=MSA_parser::check_msa_licit_and_size($inMsaFile,"fasta","no"); #HAIM ADD
	if ($ans[0] eq "OK"){$MSA_Depth=$ans[1];}    #HAIM ADD
	else {return "printColoredAlignment: ".join (" ",@ans);} #HAIM ADD
	# Print HTML start
	my %msaColors = ();
	my %msaPrintColors = ();
	my $lineCounter;
	my @line;
	my $key;
	
	my @msaRightOrder=0;
	my $msaRightOrderCounter=0;
	
	my @colorstep = (); #color steps
	$colorstep[0] = "Score1"; #Not confident
	$colorstep[1] = "Score2";
	$colorstep[2] = "Score3";
	$colorstep[3] = "Score4";
	$colorstep[4] = "Score5"; #average
	$colorstep[5] = "Score6";
	$colorstep[6] = "Score7";
	$colorstep[7] = "Score8";
	$colorstep[8] = "Score9"; #Most confident
	$colorstep[9] = "Score9"; #Most confident (the score is exactly 1)
	
	my @colorstep_code = (); #color steps
	$colorstep_code[0] = "#10C8D1"; #Not confident
	$colorstep_code[1] = "#8CFFFF";
	$colorstep_code[2] = "#D7FFFF";
	$colorstep_code[3] = "#EAFFFF";
	$colorstep_code[4] = "#FFFFFF"; #average
	$colorstep_code[5] = "#FCEDF4";
	$colorstep_code[6] = "#FAC9DE";
	$colorstep_code[7] = "#F07DAB";
		$colorstep_code[8] = "#A02560"; #Most confident
	
	# get Align max_seq_length
	my $seq = $aln->get_seq_by_pos(1);
	my $Align_width = $seq->length();
	
	# counts how many times we print the whole section (relevants to sequences longer than the sequenceLengthForDisplay) 
	for(my $blockStart=1; $blockStart<$aln->length; $blockStart+=$sequenceLengthForDisplay) {
		my $blockEnd = $blockStart+$sequenceLengthForDisplay;
		$blockEnd = $aln->length if ($blockEnd > $aln->length);
		
		# Iterate over sequences and print up to sequenceLengthForDisplay residues
		for(my $i=1;$i<=$MSA_Depth;$i++){		    #HAIM ADD
			my $seq = $aln->get_seq_by_pos($i); #HAIM ADD
			
			# NEW
			my $Color_Class="";
			# Print seq
			my @seq = split //, $seq->subseq($blockStart, $blockEnd);;
			my $gaps=0;
			for(my $pos=0; $pos<@seq; $pos++) {
				my $res = $seq[$pos];
				my $Color_Class="";
				my $prob="";
				if ($res eq '-')
				{
					$gaps++;
				}
				elsif (($res ne '-') and ($scores{$i}[$pos+1] ne "nan")) {
					$Color_Class = $colorstep[ int(9 * $scores{$i}[$pos+1]) ];
					$prob=$scores{$i}[$pos+1];
					if ($Color_Class ne "Score5")
					{
#						print JALVIEW_FEATURES "$prob\t$Code_Names{$seq->id}\t-1\t",$pos+1-$gaps,"\t",$pos+1-$gaps,"\t$Color_Class\t$prob\n";
						print JALVIEW_FEATURES "$prob\tID_NOT_SPECIFIED\t",$i-1,"\t",$pos+1-$gaps,"\t",$pos+1-$gaps,"\t$Color_Class\t$prob\n"; # COLOR JALVIEW by SEQ NUM - not SEQ ID 
					}
				}
				elsif (($res ne '-') and ($scores{$i}[$pos+1] eq "nan")) {
					$Color_Class="ScoreNaN";
#					print JALVIEW_FEATURES "$prob\t$Code_Names{$seq->id}\t-1\t",$pos+1-$gaps,"\t",$pos+1-$gaps,"\t$Color_Class\t$prob\n";
					print JALVIEW_FEATURES "NA\tID_NOT_SPECIFIED\t",$i-1,"\t",$pos+1-$gaps,"\t",$pos+1-$gaps,"\t$Color_Class\t$prob\n";      # COLOR JALVIEW by SEQ NUM - not SEQ ID 

				}
			}
		}
	}
	
	print JALVIEW_FEATURES "ENDGROUP\tGUIDANCE\n";
	close (JALVIEW_FEATURES);
}

sub make_JalView_output
{
	my $JalView_Applet_Page=shift;
	my $WorkingDir=shift;
	my $http=shift;
	# Colored MSA
	my $inMsa=shift;                      # MSA file with codes
	my $inMsa_With_names=shift;                      # MSA file with codes
	my $scores=shift;                     # hash of scores
	my $outJalviewFeaturesFile=shift;
	my $NamesCodeFile=shift;

	# Annotation Graph
	my $Jalview_AnnotFile=shift; # OUT FILE
	my $Data_File=shift;         # file ordered by the X's, first field must be X
	my $Y_label=shift;           # The Y label
	
	open (JALVIEW,">$JalView_Applet_Page") ||  return "Can't open Jalview output page: '$JalView_Applet_Page' $!";
	make_Jalview_Color_MSA("$WorkingDir$inMsa",$scores,$WorkingDir.$outJalviewFeaturesFile,$NamesCodeFile);
	make_Jalview_AnnotationGraph($WorkingDir.$Jalview_AnnotFile,$Data_File,$Y_label);
	print JALVIEW "<HTML>\n";
	print JALVIEW "<applet	CODEBASE=\"http://guidance.tau.ac.il/\"\n";
	print JALVIEW "CODE=\"jalview.bin.JalviewLite\" width=100% height=100%\n";
	print JALVIEW "ARCHIVE=\"jalviewApplet.jar\">\n";
	print JALVIEW "<param name=\"file\" value=\"$http".$inMsa_With_names."\">\n";
	print JALVIEW "<param name=\"features\" value=\"$http".$outJalviewFeaturesFile."\">\n";
	print JALVIEW "<param name=\"annotations\" value=\"$http".$Jalview_AnnotFile."\">\n";
	print JALVIEW "<param name=\"application_url\" value=\"http://www.jalview.org/services/launchApp\">\n";
	print JALVIEW "<param name=\"showbutton\" VALUE=\"false\">\n";
	print JALVIEW "<param name=\"showConservation\" VALUE=\"false\">\n";
	print JALVIEW "<param name=\"showQuality\" VALUE=\"false\">\n";
	print JALVIEW "<param name=\"showConsensus\" VALUE=\"false\">\n";
	print JALVIEW "</APPLET>\n";
	print JALVIEW "</HTML>\n";
	
}

sub validate_Seqs{
	my $workingDir=shift;
	my $in=shift;
	my $SeqType=shift;   # AminoAcids | Nucleotides | Codons
	my $MSA=shift;    # Yes,No  
	my $CodonTable=shift; # requiered if SeqType is Codons
	open (IN,$workingDir.$in) || return ('sys_error', "Validate_Seqs:Can't open '$workingDir$in': $!");
	my $seq="";
	my $seq_name="";
	my $seq_length=0;
	my $Warnning="";
	my $Counter=0;
	open (OUT,">$workingDir"."$in".".FIXED") || return ('sys_error', "Validate_Seqs:Can't open '$workingDir$in".".FIXED': $!");
	while (my $line=<IN>)
	{
		chomp ($line);
		$line=~ s/^\s+|\s+$//g;
		if (($line!~/>/)and ($line ne ""))
		{
			$seq=$seq.$line;
		}
		elsif ($line=~/^>(.*)/)
		{
			# validate prev seq
			if (($seq eq "") and ($seq_name ne ""))
			{
				return ("The sequence named '$seq_name' is missing<br>");
			}
			if (($seq ne "") and ($seq_name ne ""))
			{
				# validate seq according to type
				if ($MSA eq "Yes") # Make sure alignment length equal
				{
					$seq_length=length($seq) if ($seq_length==0); # initialize the first one
					if (length($seq)!=$seq_length)
					{
						return ("The sequences of the provided MSA are not properly aligned, For example the seq: '$seq_name' does not aligned to all others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>");
					}
					if ($SeqType eq "Codons") # Make sure that in Codon Alignment there are no stop Codons and all seq are divided by 3
					{
						my $ans=validate_seq_in_CodonAlign($seq,$seq_name,$CodonTable);
						return ($ans) if ($ans ne "OK");
					}
				}
				if ($MSA eq "No")
				{
					if ($seq=~/([-]+)$/)
					{
						$seq=~s/$1//;
						$Warnning="Gap characters (-) were removed from the end of the sequences";
					}
					if ($seq=~/[-]/)
					{
						return ("Seq: named '$seq_name' contain a gap character '-' which is illigal when sequences are submited to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>");
					}
				}
				my @ans=validate_single_seq($seq_name,$seq,$SeqType);
				if ($ans[0] eq "OK")
				{
					print OUT ">$seq_name\n";  # prev seq
					print OUT "$seq\n";         # prev seq
					$Counter++;
				}
				else
				{
					return ($ans[0]);
				}
			}
			# Start new seq
			if ($line=~/^>(.*)/)
			{
				$seq_name=$1;
				$seq_name=~ s/^\s+|\s+$//g ;
				if ($seq_name eq "")
				{
					my $Seq_Num=$Counter+1;
					return ("Seq number $Seq_Num has no sequence name; Please fix and resubmit<br>");
				}
				else
				{
#					$seq_name=$1;
					$seq="";
				}
			}
		}
	}
	# validate last seq
	if (($seq eq "") and ($seq_name ne ""))
	{
		return ("The sequence named '$seq_name' is missing<br>");
	}
	else
	{
		# validate seq according to type
		if ($MSA eq "Yes") # Make sure alignment length equal
		{
			$seq_length=length($seq) if ($seq_length==0); # initialize the first one
			if (length($seq)!=$seq_length)
			{
				return ("The sequences of the provided MSA are not properly aligned, For example the seq: '$seq_name' does not aligned to all others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>");
			}
			if ($SeqType eq "Codons") # Make sure that in Codon Alignment there are no stop Codons and all seq are divided by 3
			{
				my $ans=validate_seq_in_CodonAlign($seq,$seq_name,$CodonTable);
				return ($ans) if ($ans ne "OK");
			}
		}
		if ($MSA eq "No")
		{
			if ($seq=~/([-]+)$/)
			{
				$seq=~s/$1//;
				$Warnning="Gap characters (-) were removed from the end of the sequences";
			}
			if ($seq=~/[-]/)
			{
				return ("Seq: named '$seq_name' contain a gap character '-' which is illigal when sequences are submited to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>");
			}
		}
		my @ans=validate_single_seq($seq_name,$seq,$SeqType);
		if ($ans[0] eq "OK")
		{
			print OUT ">$seq_name\n";   # prev seq
			print OUT "$seq\n";         # prev seq
			$Counter++;
		}
		else
		{
			return ($ans[0]);
		}
	}
	#if ($Counter<4)
	#{
	#	return ("Only $Counter sequences were provided, however at least 4 sequences are requiered");
	#}
	close (OUT);
	close (IN);
	return ("OK",$Warnning,$in.".FIXED",$Counter);
}

sub validate_single_seq
{
	my $Seq_Name=shift;
	my $seq=shift;
	my $seq_type=shift;
	if (($seq!~/[ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz]+/) and ($seq_type eq "AminoAcids"))
	{
		return ("Seq: '$Seq_Name' is empty<br>");
	}
	elsif (($seq!~/[ACTGUNactgun]+/) and ($seq_type ne "AminoAcids"))
	{
		return ("Seq: '$Seq_Name' is empty<br>");
	}
	if (($seq=~/([^ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz-])/) and ($seq_type eq "AminoAcids"))#Maybe allow: _*-?
	{
		return ("Seq: '$Seq_Name' contained the character: '$1' which is not a standard Amino Acid<br>");
	}
    #----------- Amit -------------
	if (($seq=~/([^ACGTRYWSMKHBVDNUXacgtrywsmkhbvdnux-])/) and ($seq_type ne "AminoAcids"))#Maybe allow: _*-?
	{
		my $wrong_char=$1;
		if (($seq=~ /[Uu]/) and ($seq_type eq "Nucleotides"))
		{
			return ("Currently GUIDANCE does not accept 'U's in nucleotide sequences, you may consider replacing the 'U's by 'T's and re-submit. <br> In addition, seq: '$Seq_Name' contained the character: '$wrong_char' which is not a standard Nucleotide <br>");
		}
		return ("Seq: '$Seq_Name' contained the character: '$wrong_char' which is not a standard Nucleotide<br>");
	}
	if (($seq=~ /[Uu]/) and ($seq_type eq "Nucleotides"))#Maybe allow: _*-?
	{
		return ("Currently GUIDANCE does not accept 'U's in nucleotide sequences, you may consider replacing the 'U's by 'T's and re-submit.<br>");
	}
    #----------- Amit -------------
	return ("OK");
}
sub validate_seq_in_CodonAlign
# Make sure no stop Codons
{
	my $DNASequence = shift;
	my $DNASequenceName = shift;
	my $codonTableIndex = shift;
	my $stopCodon_Found="NO";
	my $AASeq="";
	my $codonTable_obj  = Bio::Tools::CodonTable -> new ( -id => $codonTableIndex );
	chomp ($DNASequence);
	my $seq_length = length($DNASequence);
    my $i =0;
	return ("Sequence '$DNASequenceName' is not a valid coding sequence: the sequence is of length $seq_length which it is not divided by 3") if ($seq_length % 3>0);
	while ($i<$seq_length-2)
	{
        $codon = substr($DNASequence, $i, 3);
        if ($codon eq '---')
		{
            $AA = '-';
        }
        else
		{
            $AA = $codonTable_obj->translate($codon);
		}
		$AASeq.= $AA;
		$i+=3;
	}
	if ($AASeq =~ m/\*/){
		return ("Sequence: '$DNASequenceName' contains a stop codon, please remove all stop codons (from all sequences) and submit to GUIDANCE again");
	}
	return ("OK");
}

1;
