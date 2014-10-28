#!/usr/bin/perl -w
# v0.2

# Runs standalone FFPred 2.0 jobs.
# The subversion (date of SVM training) must be specified by $subversion below.
# 
# The user may need to modify $FFPred_dir and/or $grid_engine_task_ID_variable
# as described below ("IMPORTANT" messages).
# 
# Also, if the temporary files produced by FFPred are not needed, they can be
# deleted automatically by removing the comment sign (#) from the appropriate
# line - look for the "Use the following line to automatically delete all
# temporary files" message below.
# 
# This script currently takes one input text file at a time (see "Sample usage"
# below), containing one or more protein sequences in FASTA format.

use strict;
use Getopt::Long;
use File::Copy;



my $subversion = '04 July 2012';

# IMPORTANT : set the following variable ($FFPred_dir) to the absolute path to
# this script.
my $FFPred_dir = '/scratch2/floris/function_prediction/FFPred2';

# IMPORTANT : if you are using this script to submit batch jobs to a
# cluster/HPC facility via e.g. SGE or other Grid Engine systems, AND you are
# submitting 'array' jobs using for example 'qsub -t' or similar commands, you
# need to set the following variable ($grid_engine_task_ID_variable) to the
# name of the environment variable that the system creates to store the index
# number of the current array job task.
# For example for the Sun Grid Engine, this environment variable is called
# 'SGE_TASK_ID', which is the default here.
# See the man page for your submission command (e.g. 'man qsub' when using SGE)
# to find out the appropriate name.
my $grid_engine_task_ID_variable = 'SGE_TASK_ID';
print "\n***Start BlastFFPred***\n";
# Sample usage :
# 
# ./FFPred.pl -i in/test.fsa -o out
# ./FFPred.pl -a -i <FFPred_main_dir>/in -o <FFPred_main_dir>/out
# 
# Check that input files/folders exist.
# Also note that if using Grid Engine systems or similar systems to submit
# batch jobs, you will need to specify full paths to all files and folders.
my $usage = "\n$0 -i [FASTA_input_file] -o [output_folder (default: \$FFPred_dir/out/)]" . 
            "\n$0 -a -i [FASTA_input_folder] -o [output_folder (default: \$FFPred_dir/out/)]\n\n";

my ($arrayjob, $fasta, $dirOut);

my $args = GetOptions(
                       "a"   => \$arrayjob,
                       "i=s" => \$fasta,
                       "o=s" => \$dirOut
                     );

die "$usage" unless (defined($fasta));

if (defined($arrayjob))
{
    die "\nBad name for the task ID variable - ABORTING !\n\n" unless (exists $ENV{$grid_engine_task_ID_variable});
    my @jobarray = sort glob("$fasta/*");
    $fasta = $jobarray[$ENV{$grid_engine_task_ID_variable}-1];
}

$dirOut = "$FFPred_dir/out" unless (defined($dirOut));
mkdir $dirOut unless (-d $dirOut);
die "$usage" unless ((-T $fasta) && (-w $dirOut));

# Initialise useful variables.
my $submit_datetime = Timestamp();

my $jobs_directory = ReadConfigPATH($FFPred_dir)
    or die "\n----- error in sub ReadConfigPATH -----\n\n";

my $GOdomains = {
                  'biological_process' => 'BP',
                  'molecular_function' => 'MF'
                };

# Read input sequences into a hashref.
my $inputs = {};
ReadFasta($fasta, $inputs)
    or die "\n----- error in sub ReadFasta -----\n\n";

# Perform FFPred jobs.
foreach my $id (keys %$inputs)
{
    my ($error, $job, $pred) = (0, {}, {});
    
    $job->{'id'} = $id;
    $job->{'seq'} = $inputs->{$id};
    $job->{'len'} = length($job->{'seq'});
    $job->{'submitted'} = $submit_datetime;
    
    # Find a root name for output files using this sequence's FASTA header.
    my $output_name = my $output_rootname = ($job->{'id'} =~ /([\w\-]{1,10})/) ? $1 : 'default';
    my $output_folder = "$dirOut/FFPred_${output_rootname}";
    my $suffix_dir = 0;
    while (!mkdir($output_folder))
    {
        die "Too many identical output folder names !\n" if (++$suffix_dir > 100000);
        $output_name = "${output_rootname}_${suffix_dir}";
        $output_folder = "$dirOut/FFPred_${output_name}";
    }
    
    $job->{'out'} = "$output_folder/$output_name";
    
    GetMD5($job);
    
    print STDERR "\n".Timestamp()." - Started job $job->{'id'}.\n" . 
                 "Input file $fasta\n" . 
                 "Input sequence's md5 code $job->{'md5'}\n" . 
                 "Output folder $output_folder\n";
    
    if ($job->{'len'} < 15)
    {
        print STDERR "\n\nBEWARE !\n" . 
                     "$job->{'id'} was NOT processed !!!\n" . 
                     "Its sequence is too short: minimum 15aa, $job->{'len'} aa submitted.\n\n\n";
        next;
    }
    elsif ($job->{'seq'} =~ /^[CTAG]+$/)
    {
        print STDERR "\n\nBEWARE !\n" . 
                     "$job->{'id'} was processed, but it looks like DNA !!!\n" . 
                     "FFPred is intended for proteins.\n\n\n";
    }
    
    MkDir($jobs_directory, $job);
    
    print STDERR "\n".Timestamp()." - Temporary folder $job->{'dir'} - Started featurama.\n\n";
    $error += RunFeaturama($FFPred_dir, $job);						# Run features.pl in the featurama dir which is located in FFPred dir
    print STDERR "\n".Timestamp()." - Temporary folder $job->{'dir'} - Finished featurama!\n";
    
    die "Error ($error) while running featurama - ABORTING !\n" if ($error);   
     
    # Use the following line to automatically delete all temporary files.
    #Clean_temp($job);
    
    print STDERR "\n".Timestamp()." - Finished job $job->{'id'}.\n\n\n";
}



################################################## Subroutines after this line.



sub ReadConfigPATH
{
    my ($rootdir) = @_;
    my $cfgPATH = "";
    
    open(CONFIG, "<", "$rootdir/CONFIGBLAST") or print STDERR "Cannot read CONFIGBLAST - $!\n" and return;
    
    while (defined(my $line = <CONFIG>))
    {
        next if ($line =~ /^#/);
        if ($line =~ /\bPATH\s+(\S+)/)
        {
            $cfgPATH = $1;
            $cfgPATH = "$rootdir/$cfgPATH" unless (substr($cfgPATH, 0, 1) eq '/');
        }
    }
    
    close(CONFIG) or print STDERR "Cannot close CONFIGBLAST - $!\n";
    
    print STDERR "'PATH' for temporary jobs not found in CONFIGBLAST file!\n" unless ($cfgPATH);
    return $cfgPATH;
}

sub ReadFasta
{
    my ($fasta_filepath, $inputs_hash) = @_;
    
    local $/;
    open(FASTA, "<", $fasta_filepath) or print STDERR "Cannot open file $fasta_filepath - $!\n" and return;
    my $content = <FASTA>;
    close(FASTA) or print STDERR "Cannot close file $fasta_filepath - $!\n";
    
    print STDERR "Invalid input FASTA file !\n" and return unless ($content =~ />/);
    $content =~ s/\A.*?^>/>/ms; # Avoid any initial lines without a valid FASTA header.
    my $key;
    
    foreach my $line (split("\n", $content))
    {
        if ($line =~ /^>/)
        {
            $line =~ s/^>\s*//;
            $line =~ s/\s*$//;
            
            if ($line ne "")
            {
                $key = $line;
            }
            else
            {
                print STDERR "Invalid FASTA header in input file !\n" and return;
            }
        }
        elsif ($line !~ /^\s*$/)
        {
            $line =~ s/\s//g;
            $line = uc $line;
            
            if ($line =~ /^[A-Z]+$/)
            {
                $inputs_hash->{$key} .= $line;
            }
            else
            {
                print STDERR "Invalid line in input file !\n" and return;
            }
        }
    }
    
    return 1;
}

sub GetMD5
{
    my ($job_hash) = @_;
    
    use Digest::MD5 qw(md5_hex);
    
    my %hextobinhash =   ('x0'=>'0000', 'x1'=>'0001', 'x2'=>'0010',
                          'x3'=>'0011', 'x4'=>'0100', 'x5'=>'0101',
                          'x6'=>'0110', 'x7'=>'0111', 'x8'=>'1000',
                          'x9'=>'1001', 'xa'=>'1010', 'xb'=>'1011',
                          'xc'=>'1100', 'xd'=>'1101', 'xe'=>'1110',
                          'xf'=>'1111');
    
    my %bintohex32hash = ('x00000'=>'0', 'x00001'=>'1', 'x00010'=>'2',
                          'x00011'=>'3', 'x00100'=>'4', 'x00101'=>'5',
                          'x00110'=>'6', 'x00111'=>'7', 'x01000'=>'8',
                          'x01001'=>'9', 'x01010'=>'a', 'x01011'=>'b',
                          'x01100'=>'c', 'x01101'=>'d', 'x01110'=>'e',
                          'x01111'=>'f', 'x10000'=>'g', 'x10001'=>'h',
                          'x10010'=>'i', 'x10011'=>'j', 'x10100'=>'k',
                          'x10101'=>'l', 'x10110'=>'m', 'x10111'=>'n',
                          'x11000'=>'o', 'x11001'=>'p', 'x11010'=>'q',
                          'x11011'=>'r', 'x11100'=>'s', 'x11101'=>'t',
                          'x11110'=>'u', 'x11111'=>'v');
    
    my $enc = md5_hex($job_hash->{'seq'});
    my $md5cut = substr($enc, -20, 20);
    
    my ($char, $binary, $hex_md5, $bin);
    
    for (my $x = 0; $x < 20; $x++)
    {
        $char = "x" . substr($md5cut, $x, 1);
        $binary .= $hextobinhash{$char};
    }
    
    $hex_md5 = "";
    
    for (my $y = 0; $y < 80; $y += 5)
    {
        $bin = "x" . substr($binary, $y, 5);
        $hex_md5 .= $bintohex32hash{$bin};
    }
    
    $job_hash->{'md5'} = $hex_md5;
}

sub MkDir
{
    my ($jobs_dir, $job_hash)= @_;
    
    $job_hash->{'dir'} = my $job_root = "$jobs_dir/" . substr($job_hash->{'md5'}, 0, 5);
    
    # This avoids race conditions when many jobs are run at the same time.
    my $suffix_tempdir = 0;
    while (!mkdir($job_hash->{'dir'}))
    {
        die "Too many identical temporary folder names !\n" if (++$suffix_tempdir > 100000);
        $job_hash->{'dir'} = "${job_root}_${suffix_tempdir}";
    }
    
    open(OUT, ">", "$job_hash->{'dir'}/$job_hash->{'md5'}.fsa");
    print OUT ">$job_hash->{'md5'}\n$job_hash->{'seq'}\n";
    close(OUT);
}

sub Timestamp
{
    my @MONTHS = ('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');
    my ($sec, $min, $hour, $day, $month, $year) = (localtime)[0..5];
    
    return "$day $MONTHS[$month] " . ($year+1900) . ", $hour:$min:$sec";
}

sub RunFeaturama
{
    my ($rootdir, $job_hash) = @_;
    print STDERR "**Call Blast_features.pl** \n";
    my $featurama = "$rootdir/featurama/Blast_features.pl " .
                    "-d $job_hash->{'dir'} " . 
                    "-i $job_hash->{'dir'}/$job_hash->{'md5'}.fsa " . 
                    "-o $job_hash->{'dir'}/$job_hash->{'md5'}.results " . 
                    "-fconfig $job_hash->{'dir'}/$job_hash->{'md5'}.featcfg";
    $featurama .= " -orphan" if ($job_hash->{'len'} <= 20);
    
    return system($featurama);
}

sub Clean_temp
{
    my ($job_hash) = @_;
    
    unlink(glob("$job_hash->{'dir'}/$job_hash->{'md5'}" . "[._]*"));
    rmdir $job_hash->{'dir'} or print STDERR "Could not delete folder $job_hash->{'dir'} - $!\n";
}
