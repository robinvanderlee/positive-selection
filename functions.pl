#!/usr/bin/env perl
#####################################
## Robin van der Lee               ##
## robinvanderlee AT gmail DOT com ##
############################################################################################################
## Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts ##
## Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen                                ##
############################################################################################################
use strict;
use warnings;

### PRINTING
sub newline{
	if(defined $_[0]){
		print "\n" x $_[0];
	} else {
		print "\n";
	}
}

sub nl{
	newline(@_);
}

sub printArray{
	foreach(@_){
		print;
		newline();
	}
}

sub printArray2{
	foreach(@_){
		print $_ . "\t";
	}
	print "\n";
}

sub printArrayofHashes{
	foreach(@_){
		my %tmphash = $_;
		printHash(\%tmphash);
		separator();
	}
	print "\n";
}

sub printHash{
	my %hash = %{ $_[0] };
	foreach(keys %hash){
		print $_ . " - " . $hash{$_} . "\n";
	}
}

sub printHashOfArrays{
	my %hash = %{ $_[0] };
	foreach(keys %hash){
		my @array = @{ $hash{$_} };
		print $_ . " - [ " . join(",", @array) . " ]\n";
	}
}

### SEPARATORS
sub separator{
	print generic_separator("=");
}

sub separator2{
	print generic_separator("-");
}

sub separator3{
	print generic_separator("~");
}

sub separator4{
	print generic_separator("+");
}

sub separator5{
	print generic_separator("\"");
}

sub separator6{
	print generic_separator("#");
}

sub generic_separator{
	my $character = shift;
	my $separator_string = $character x 60 . "\n";
	return $separator_string;
}

### TIME
sub getTimeAndDate{
	my @Month_name = ( "January", "February", "March", "April", "May", "June",
  "July", "August", "September", "October", "November", "December" );
	( my $sec, my $min, my $hour, my $day, my $month, my $year ) = ( localtime ) [ 0, 1, 2, 3, 4, 5 ];
	printf "%02d:%02d:%02d %02d %s %04d\n", $hour, $min, $sec, $day, $Month_name[$month], $year+1900;
}

# sub getTimestamp{
#    ( my $sec, my $min, my $hour, my $day, my $month, my $year ) = ( localtime ) [ 0, 1, 2, 3, 4, 5 ];
# 	return sprintf "%04d%02d%02d-%02d%02d", $year+1900, $month + 1, $day, $hour, $min;
# }

sub getTimestampSec{
   ( my $sec, my $min, my $hour, my $day, my $month, my $year ) = ( localtime ) [ 0, 1, 2, 3, 4, 5 ];
	return sprintf "%04d%02d%02d-%02d%02d%02d", $year+1900, $month + 1, $day, $hour, $min, $sec;
}

sub getTimestampMin{
   ( my $sec, my $min, my $hour, my $day, my $month, my $year ) = ( localtime ) [ 0, 1, 2, 3, 4, 5 ];
	return sprintf "%04d%02d%02d-%02d%02d", $year+1900, $month + 1, $day, $hour, $min;
}

sub getTimestampHour{
   ( my $sec, my $min, my $hour, my $day, my $month, my $year ) = ( localtime ) [ 0, 1, 2, 3, 4, 5 ];
	return sprintf "%04d%02d%02d-%02d", $year+1900, $month + 1, $day, $hour;
}

sub getTimestampDay{
   ( my $sec, my $min, my $hour, my $day, my $month, my $year ) = ( localtime ) [ 0, 1, 2, 3, 4, 5 ];
	return sprintf "%04d%02d%02d", $year+1900, $month + 1, $day;
}

return 1;
