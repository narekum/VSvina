#!/usr/bin/perl 

################################################################################
# Parse sdf file from drugbank
# Written by Narendra Kumar
# 15th April 2020
#
# narekum@gmail.com
# https://narekum.github.io
#############################

BEGIN {
        use Cwd qw(realpath cwd);
        use File::Basename;
        our ($fn, $dir) = fileparse(realpath($0));
        $dir =~ s/\/bin\/$// ;
}


################################################################################
# Import modules
use File::Path qw(make_path);
use LWP::Simple;
use Getopt::Long ;

print STDERR '
__     ______                            
\ \   / / ___| _ __   __ _ _ __ ___  ___ 
 \ \ / /\___ \| "_ \ / _" | "__/ __|/ _ \
  \ V /  ___) | |_) | (_| | |  \__ \  __/
   \_/  |____/| .__/ \__._|_|  |___/\___|
              |_|                        
    SDF parsing for Virtual Screening   
';


my $help ;
my $sdf ;
my $fields ;
#my $report ;

################################################################################
# Get commandline arguments

GetOptions ("h|help"=>\$help,                     # --help         : print this help
            "s|sdf=s"=>\$sdf,                     # --sdf          : sdf input file
            "f|fields=s"=>\$fields,               # --fields       : comma separated fields to report in output
            "r|report!" =>\$report                # --report       : reports the fields in the input sdf file 
           )
or &PRINTHELP($fn);

################################################################################
# Check input and warn

if( defined $help || !$sdf ) {
        &PRINTHELP($fn);
        exit;
}


################################################################################
# Report the available fields in the given sdf file

if (defined $report ) {

	print "\n";
	print "    Input sdf file \"$sdf\" has the following fields:\n\n";

	my $all_fields=REPORT_FIELDS($sdf);
	print "    $_\n" foreach sort keys %$all_fields ;

	print "\n\n\n";

	exit;
}

################################################################################
# Report user define fields in the output. if -f option is true

my @default_fields = ("DRUGBANK_ID","GENERIC_NAME");

my @fields=();
if ( defined $fields ) {
	@fields = split ( /,|;|:/ , $fields ) ;
}

# Remove from user supplied fields the default fields.
my @fields_to_print=();
my @extra_fields=();
my $element ;
OUTER: foreach $element (@fields) {
	INNER: foreach (@default_fields) {
		if ( $element eq $_ ) {
			next OUTER;
		}
	}
	push @extra_fields,$element ;
} 

# Field to print in that order
#@fields_to_print = (@default_fields,@extra_fields);
push @default_fields,@extra_fields ;
print "$_\n" foreach @fields_to_print ;

my @prefix_ext = ("structures",".pdbqt");

$/ = '$$$$';

my $count=0;

open my $fh, $sdf ;

while  (<$fh>) {
	next if $_ =~ /^$/ ;
	$count++;
	my $report_table =PRINT_FILEDS($count,$_,\@default_fields,\@prefix_ext);
	print "$report_table\n";
}

close $fh ;

$/="\n";

################################################################################
# Report the given field values in that order

#my @fields=("DATABASE_ID", "INCHI_KEY", "FORMULA");
sub PRINT_FILEDS {
	my ($record_number, $record, $fields_ref, $prefix_ext_ref)=@_;
	my @fields =  @$fields_ref ;
	my @prefix_ext = @$prefix_ext_ref ;
	my @record_lines = split(/\n/, $record) ;
	my $line_count=0;
	my %fields_values=();
	my @fields_values;
	my $report_table ;
	foreach my $lines (@record_lines) {
		chomp;
		$line_count++;
		foreach my $field_lines (@fields) {
			if ($lines =~ /$field_lines/) {
				$fields_values{$field_lines} = $record_lines[$line_count] ;
			}
		}
	}
	$record_number = join "", $prefix_ext[0] , $record_number , $prefix_ext[1];
	foreach (@fields) {
		push @fields_values,$fields_values{$_} ;
	}
	$report_table = join " ",$record_number,  @fields_values ;
	return ($report_table);
}




################################################################################
# Obtain all the extractable fields from drugbank sdf file

sub REPORT_FIELDS { 
	my ($drugbank_sdf) = @_ ;
	$/ = '$$$$' ;
	open my $fh, $drugbank_sdf ;
	my @each_record=();
	my $count=0;
	my %all_fields=();
	while (<$fh>) {
		$count++;
		@each_record=();
		@each_record = split (/\n/, $_ ) ;
		foreach (  @each_record ) {
			if ($_ =~ /^\s*>\s*<(.*?)>/ ) {  # "> <DATABASE_ID>" this is how the record looks like
				my $field_name = $1 ;
				#print "$count $field_name\n";
				$all_fields{$field_name}=1;
			}
		}
	}
	close $fh;
	$/ = '\n' ;
	return (\%all_fields);
}


################################################################################
# Pring help routine 

sub PRINTHELP {
        my ($fn)=@_;

        print <<DES;

Usage: $fn -s drugbank.sdf 

Options:

  -h, --help      Print this help
  -s, --sdf       Sdf input file.
  -f, --fields    comma separated list of fields to report in output
  -r, --report    reports the fields in the input sdf file

Examples:
\$ perl $fn -s drugbank.sdf -r  
\$ perl $fn -s drugbank.sdf -report  

\$ perl $fn -s drugbank.sdf 
\$ perl $fn -s drugbank.sdf -f EXACT_MASS,INCHI_KEY

DES

exit;
}
