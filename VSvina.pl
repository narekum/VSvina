#!/usr/bin/perl 

################################################################################
# Screening with Autodock Vina
# Written by Narendra Kumar
# 14th April 2020
#
# narekum@gmail.com
# https://narekum.github.io
#############################


BEGIN {
        use Cwd qw(realpath cwd);
        use File::Basename;
        our ($fn, $dir) = fileparse(realpath($0));
}

#use strict ;
#use LWP::Simple;
use Getopt::Long;
use Cwd qw(cwd);
use File::Copy;

print '
__     ______       _             
\ \   / / ___|_   _(_)_ __   __ _ 
 \ \ / /\___ \ \ / / | \'_ \ / _\' |
  \ V /  ___) \ V /| | | | | (_| |
   \_/  |____/ \_/ |_|_| |_|\__._|
    Virtual Screening of Ligands                                  

' ;


$vina="~/local/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina";

GetOptions ('h|help'=>\$help,                     # --help       : print this help
            "l|liganddir=s" => \$ldir,            # --liganddir  : directory containing pdbqt files of the ligands 
            "r|receptor=s" => \$rpdbqt,           # --receptor   : receptor pdbqt file
            "o|outdir=s" => \$outdir,             # --outdir     : output directory (will be created)
            "c|conf=s" => \$conf,                 # --conf       : autodock vina configuration file 
            "m|mapfile=s" => \$map                # --mapfile    : mapfile maping names of ligands to its pdbqt file
           )
or &PRINTHELP($fn);

if(defined $help || !$ldir || !$rpdbqt || !$outdir || !$conf) {
        &PRINTHELP($fn);
        exit;
}

################################################################################
# Checking for the existence of --liganddir

if ( -d $ldir ) {
        print "\n   Ligand Directory found: \"$ldir\"  \n\n";
} else {
	print "\n   Lignad Directory NOT found:$ldir ";
	exit;
}

opendir my $dh, $ldir or die "\n    Could not open '$ldir' for reading: $!\n";
my @things = readdir $dh ;
closedir $dh ;

my @ldir_pdbqts = () ;
foreach my $thing (@things) {
    if ( $thing =~ /\.pdbqt$/) {
	push @ldir_pdbqts,$thing ;
    }
}

print "\n   No of ligands found in $ldir for screening: ", scalar(@ldir_pdbqts),  "\n";


################################################################################
# Checking receptor

if (-e $rpdbqt) {
	print "\n   Receptor file $rpdbqt found !!!\n\n" ;
} else {
	print "\n   Receptor file $rpdbqt NOT found !!! Exiting!\n\n" ;
	exit;
	
}


################################################################################
# Checking configuration file

if (-e $conf ) {
	print "\n   Vina Configuration file $conf found !!!\n\n" ;
} else {
	print "\n   Vina Configuration file $conf NOT found !!! Exiting\n\n" ;
	exit;

}


################################################################################
# Checking mapfile 

if (-e $map ) {
	print "\n   Map file \"$map\" found !!!\n\n" ;
} else {
	print "\n   Map file \"$map\" NOT found !!! Exiting\n\n" ;
	exit;

}


################################################################################
# Checking outdirectory 

my $pwd = cwd ; 
if ( $outdir !~ /^\s*\// ) {
	$outdir =~ s/.\/// ;
	$outdir = "$pwd/$outdir";
} 

if (-d $outdir ) {
	print "\n   Output directory \"$outdir\" already exists ! Not creating it!!! \n\n" ;
	#exit;
} else {
	print "\n   Creating output directory \"$outdir\" \n\n" ;
	mkdir $outdir ;

}


################################################################################
# Creating results, log directories

if (! -d "$outdir/logs" ) { mkdir "$outdir/logs" ; } else { print "        \"$outdir/logs\" already exists! Not creating it!!!\n" }  
if (! -d "$outdir/results" ) { mkdir "$outdir/results" ; } else { print "        \"$outdir/results\" already exists! Not creating it!!!\n" }  

################################################################################
# Copying files to output directory

copy($rpdbqt, $outdir ) or die ("   Copying receptor file \"$rpdbqt\" to output directory failed: $!") ;
copy($map, $outdir ) or die ("   Copying map file \"$map\" to output directory failed: $!") ;
copy($conf, $outdir ) or die ("   Copying configuration file \"$conf\" to output directory failed: $!") ;

$rpdbqt =~ s/.*\///; # remove path from before the receptor file so that the copied one may be used
$conf =~ s/.*\///;   # remove path from before the conf file so that the copied one may be used

################################################################################
# Running autodock vina

foreach my $ligand ( @ldir_pdbqts) {
	$ligand =~ s/.pdbqt\s*// ;
	print "    Preparing to  dock $ligand to $rpdbqt \n\n";
	if (! -e "$outdir/logs/$ligand.log" ) {
		print "    Running \"$vina --receptor $outdir/$rpdbqt --ligand $ldir/$ligand.pdbqt --config $outdir/conf.txt --out $outdir/results/$ligand.pdbqt --log $outdir/logs/$ligand.log\"" ;
		system ("$vina --receptor $outdir/$rpdbqt --ligand $ldir/$ligand.pdbqt --config $outdir/conf.txt --out $outdir/results/$ligand.pdbqt --log $outdir/logs/$ligand.log") ;
		print "    Finished docking $ligand to $rpdbqt \n\n\n\n" ;
	} else {
		print "    $outdir/logs/$ligand.log already exists!!! Skipping docking for \"$ligand.log\" !!! \n\n\n\n";
	}
}

################################################################################
# Compiling results

# Reading the mapfile 

open $fh, "$outdir/$map" or die ("   Can not open $mapfile for reading: $! \n");
chomp( my @map_ids=<$fh>) ;
close $fh;

my %mapped_ids = ();
my %mapped_ids_rest = ();
foreach (@map_ids) {
	my @data = split(/[\s\t]+/,$_);
	my $ligand_file_name = $data[0] ;
	my $drugbank_id = $data[1] ;
	my $drugbank_rest_data = join "\t", @data[2 .. $#data] ;
	$mapped_ids{$ligand_file_name} = $drugbank_id ;
	$mapped_ids_rest{$ligand_file_name} = $drugbank_rest_data ;
}

###############################################################################
# Reading the results

my @best_energies=();
my @header_energies=();
foreach my $ligand ( @ldir_pdbqts) {
	$ligand =~ s/.pdbqt\s*// ;
	my @temp_storage=();
	open $fh, "$outdir/logs/$ligand.log" ;
	chomp(@temp_storage=<$fh>);
	close $fh;
	my $best_energy=$temp_storage[25] ;
	my $drugbank_id = $mapped_ids{"$ligand.pdbqt"} ;
	my $drugbank_rest = $mapped_ids_rest{"$ligand.pdbqt"} ;
	@header_energies = @temp_storage[22,23,24] ;
	push @best_energies, "$drugbank_id\t$best_energy\t$drugbank_rest", ;
}
	
################################################################################
# Sorting the compounds with their energies

my @best_compounds = sort { (split /\s+/, $a )[2] <=> (split /\s+/ , $b )[2] } @best_energies ;

################################################################################
# Printing results 

open $fh, ">$outdir/${rpdbqt}_docking_results.txt" or die ("    Cant not open \"$outdir/${rpdbqt}_docking_results.txt\"for writing: $! \n\n ") ;

print $fh "       \t$_\n" foreach @header_energies ; 
print $fh "$_\n" foreach @best_compounds ;

print "\n\n   DONE SCREENING\n\n";

################################################################################

sub PRINTHELP {
        my ($fn)=@_;

        print <<DES;

Usage: $fn < -l <LIGAND_DIR> -r <RECEPTOR_PDBQT_FILE> -o <OUTPUT_DIR> -c <VINA_CONFIG_FILE>

Options:
  -h, --help       Print this help.
  -l, --liganddir  Directory containing pdbqt files of the ligands
  -r, --receptor   Receptor pdbqt file
  -o, --outdir     Output directory (will be created)
  -c, --conf       Autodock vina configuration file
  -m, --mapfile    Mapfile matching names of ligands to its pdbqt file
                   Mapfile is a tab delimited in the following format

                   NAME_OF_THE_LIGAND<TAB>FILE_NAME.pdbqt

Examples:
\$ perl $fn --liganddir ./ligands --receptor receptor.pdbqt --outdir screeing_results --conf configuration.txt 
\\or\\
\$ perl $fn -l ./ligands -r receptor.pdbqt -o screeing_results -c configuration.txt 

DES
exit;
}
