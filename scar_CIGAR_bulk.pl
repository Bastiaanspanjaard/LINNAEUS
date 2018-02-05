#!/usr/bin/perl -w

# Description
# Extract scars from paired-end fastqs in bulk sequencing data.
# The protocol generates reads where R1 has a 8-base barcode, R2 the scar. 
# Before read extraction, we take the UMI and barcode-information from the R1 and append
# it to the R2, and map the R2. The R2 is reverse mapped.

# Input and usage:
# ./scar_CIGAR_bulk.pl  -R1 [input R1]
#			-R2 [input R2]
#			-op [output prefix]
#			-t  [threads for mapping]
#			-r  [reference file to map to]
#			-bc [barcode list]
#                 	-g  [requested gene]
#			-k  [amount of bases in R2 to consider (starting from the primer)]
#                       -l  [allowed length of R2]
#                       -ps [enforced primer sequence]

# Output
# A fastq-file named [prefix].fastq that combines the UMI and barcode-information from 
# R1 (in the readname field) together with the R2 sequence and quality.
# A samfile named [prefix].sam that results from mapping this fastq to the reference.
# A file named [prefix]_scars.txt that contains all scars found in the samfile, one line
# per scar.
# A count table named [prefix]_cscars_coutc.txt in which all scar reads are assigned
# to their barcodes - this is before doing a UMI-correction.
# A count table named [prefix]_cscars_couta.txt in which all scar abundances are assigned
# to their barcodes - this is after doing a UMI-correction and this file should be used
# for all downstream analysis.
# A statistics file named [prefix]_stats.txt with the total amount of reads, the number
# of reads mapped, and the number of reads over applied thresholds (thresholds are
# currently hardcoded but can be changed in the code below).

# Created by: B. Spanjaard

# Dependencies
use strict;
use warnings;
use Getopt::Long;
use POSIX;

# Globals - could change this to input if necessary
my $barcodelength = 8;
my $barcodestart = 0;
my $min_CIGAR_reads = 20;
my $min_CIGAR_fraction = 0.05;
my $protocol = 4;

my $include_barcode_rewrite = 1;
my $include_map = 1;
my $include_barcode_read = 1;
my $include_write_scars = 1;


# Subs
sub FindEndLocation
{
    # Given the position of the leftmost matching base and the cigar score,
    # returns the genomic location of the rightmost base of the read.
    
    my $position = $_[0];
    my $cigar = $_[1];

    my @cigar_numbers = split(/[A-Z]/, $cigar);
    # Calculate the total number of bases in the cigar. This includes deleted
    # bases. Even though these are not in the read, they are in the genome
    # and should therefore be incorporated in the genomic position.
    my $frag_length_w_D = 0;
    foreach my $cigar_number (@cigar_numbers)
    {
	$frag_length_w_D += $cigar_number;
    }
    
    my @cigar_chars = split(/\d+/, $cigar);
    # The number of mismatches before the leftmost matching base should not be counted
    # in the length of the whole read. The number of mismatches after the rightmost
    # matching base should also be subtracted from the read.
    my $s_count = 0;
#	if cigar starts with S, add number belonging to first S; then, if cigar ends with S,
#	add number belonging to last S.
#	Pseudo: while $cigar =~ /([0-9]*S)/g{@x = split("S", $1); $s_count = $s_count + $x[0];}

    if($cigar_chars[1] eq "M")
    {
	$s_count = 0;
    } elsif ($cigar_chars[1] eq "S")
    {
	my @cigar_S_split = split("S", $cigar);
	$s_count = $cigar_S_split[0];
    } else
    {
	print("Unexpected cigar! First character $cigar_chars[1].\n");
    }
    my $start_location = $position - $s_count;

    my @i_split = split("I", $cigar);
    pop(@i_split);
    # The number of insertions does not have influence on the genomic end
    # location, but they do take up bases in the read.
    my $total_i_count = 0;
    foreach my $part (@i_split)
    {
	my ($i_count) = $part =~ /(\d+)$/;
    # print(" $i_count");
	$total_i_count += $i_count;
    }
    
    my $end_location = $start_location + $frag_length_w_D - 1 - 
	$total_i_count;
    
    return $end_location;
}

sub GetBarcode
{
	my $R1_bc_UMI = $_[0];
	
	my $barcode = "";
	$barcode = substr($R1_bc_UMI, 5, $barcodelength);
	
	return $barcode;
}

# Get input parameters and open files
my $in_R1_name = "";
my $in_R2_name = "";
my $out_prefix = "";
my $threads = "";
my $reference = "";
my $barcode_name = "";
my $gene = "";
my $R2_keep = "";
my $length = "";
my $primer_sequence = "";
GetOptions
(
    'R1=s' => \$in_R1_name,
    'R2=s' => \$in_R2_name,
    'op=s' => \$out_prefix,
    't=s' => \$threads,
    'r=s' => \$reference,
    'bc=s' => \$barcode_name,
    'g=s' => \$gene,
    'k=i' => \$R2_keep,
    'l=i' => \$length,
    'ps=s' => \$primer_sequence
);
open(my $R1_file, $in_R1_name)
	or die "Cannot find $in_R1_name\n";
open(my $R2_file, $in_R2_name)
	or die "Cannot find $in_R2_name\n";
open(my $fastq_out_file, ">".$out_prefix.".fastq");
open(my $barcodefile, $barcode_name)
    or die "Cannot find $barcode_name\n";

# Write R2 information and R1 barcode and UMI to new fastq
if ($include_barcode_rewrite == 1)
{
while(my $R1_header = <$R1_file>)
{
	chomp($R1_header);
	my $R1_sequence = <$R1_file>;
	chomp($R1_sequence);
	my $R1_plus = <$R1_file>;
	my $R1_qual = <$R1_file>;
	
	my $R1_UMI_bc = substr($R1_sequence, $barcodestart, $barcodelength);
	
	my $R2_header = <$R2_file>;
	chomp($R2_header);
	my @R2_header_fields = split(" ", $R2_header);
	my $R2_out_header = $R2_header_fields[0]." BC:Z:".$R1_UMI_bc;
	print $fastq_out_file "$R2_out_header\n";
	
	my $R2_sequence = <$R2_file>;
	chomp($R2_sequence);
	my $R2_out = substr($R2_sequence, 0, $R2_keep);
	print $fastq_out_file "$R2_out\n";
	my $R2_plus = <$R2_file>;
	chomp($R2_plus);
	print $fastq_out_file "$R2_plus\n";
	my $R2_qual = <$R2_file>;
	chomp($R2_qual);
	my $R2_qual_out = substr($R2_qual, 0, $R2_keep);
	print $fastq_out_file "$R2_qual_out\n";
}
}
close $fastq_out_file;

# Map new R2
my $sam_name = $out_prefix.".sam";
if ($include_map == 1)
{
my $map_command = "bwa mem -t ".$threads." -v 0 -C ".$reference." ".$out_prefix.".fastq > ".$sam_name;
system($map_command) == 0 or die "Could not execute ".$map_command."\n";
}

# Read samfile and write CIGAR-scars to output-file
my %CIGAR_hash = (); # Counts all reads in a CIGAR including duplicates, keys are CIGARs. 
my %scar_UMI_hash = (); # Counts all reads including duplicates, keys are 
# CIGAR_sequence_location, barcode and UMI.
my %scar_hash = (); # Counts all reads including duplicates, keys are CIGAR_sequence_location and
# barcode.
my %scar_UMI_unique_hash = (); # Counts all unique reads, keys are
# CIGAR_sequence_location and barcode.
my $reads = 0;
my $mapped_reads = 0;
my $mapped_reads_barcode = 0;
my $mapped_reads_length = 0;
my $mapped_reads_with_PCRp = 0;
my %primer_hash = ($primer_sequence => 1,);
my $primer_length = length($primer_sequence);
my $mapped_reads_length_with_PCRp =0;
if ($include_write_scars == 1)
{
my $scarfile_out_name = $out_prefix."_scars.txt";
open (my $scarfile_out, ">".$scarfile_out_name);
open(my $samfile_in, $sam_name)
    or die "Cannot find $sam_name\n";
print $scarfile_out "CIGAR\tBarcode\tLeft location\tFlag\tSNP-code\tSequence\n";
while (my $samline = <$samfile_in>)
{
    # Find read
    chomp($samline);
    if(substr($samline, 0 , 1) eq "\@")
    {
	next;
    }
    my @R1_fields = split(" ", $samline);
    if($R1_fields[1] > 256)
    {
	next;
    }
    $reads++;

    # Search for the barcode-UMI comment.
    my $bc_UMI_index = 13;
    while(substr($R1_fields[$bc_UMI_index], 0, 5) ne "BC:Z:")
    {
    	$bc_UMI_index++;
    }
    
    my $barcode = GetBarcode($R1_fields[$bc_UMI_index]);
    
    # Filter mapped reads
    if(!($R1_fields[2] eq $gene))
    {
    	next;
    }
    $mapped_reads++;
    
    # Filter reads with correct length.
    my $read_length = "";
   	$read_length = length($R1_fields[9]);
    if($read_length != $length)
    {
    	next;
    }
    $mapped_reads_length++;

    # Filter reads with correct primer.
    my $read_primer = "";
    $read_primer = substr($R1_fields[9], -$primer_length);
    if(!exists($primer_hash{$read_primer}))
    {
    	next;
    }
    $mapped_reads_length_with_PCRp++;

    # Add scar to tables.
    my $cigar = $R1_fields[5];
    my $location = $R1_fields[3];
    my $flag = $R1_fields[1];
    my $SNP = $R1_fields[12];
    my $sequence = $R1_fields[9];
    	
    print $scarfile_out "$cigar\t$barcode\t$location\t$flag\t$SNP\t$sequence\n";
}
close($samfile_in);
close($scarfile_out);
}

# Write statistics file
my $statfile_out_name = $out_prefix."_stats.txt";
open (my $stat_outfile, ">".$statfile_out_name);
print $stat_outfile "Reads: $reads\n";
my $reads_map_perc = 100 * $mapped_reads/$reads;
print $stat_outfile "Mapped to $gene: $mapped_reads ($reads_map_perc %)\n";
my $reads_map_length_perc = 100 * $mapped_reads_length/$reads;
print $stat_outfile "Mapped to $gene with right length: $mapped_reads_length ($reads_map_length_perc %)\n";
my $reads_map_length_PCRp_perc = 100 * $mapped_reads_length_with_PCRp/$reads;
print $stat_outfile "Mapped to $gene with right length and primer: $mapped_reads_length_with_PCRp ($reads_map_length_PCRp_perc %)\n";close($stat_outfile);
