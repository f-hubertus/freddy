#!/usr/bin/perl -w

package Seq2Mtx;

use strict;
use File::Copy;
use base 'Tool';

print STDERR "***Start Seq2MTX-Tool*** \n";

sub new
{
    print STDERR "**Seq2 new** \n";
    my ($class, $exe_path, $pfilt_path, $prefix, $fasta_path) = @_;
    
    my $self = {
                 'exe'        => "$exe_path/seq2mtx",
                 'pfilt'      => "$pfilt_path/pfilt",
                 'infile'     => "$prefix.fsa",
                 'outfile'    => "$prefix.mtx",
                 'fasta_path' => $fasta_path,
                 'err'        => 0
               };
    
    bless $self, $class;
    
    $self->{'pfilt_cmd'} = "$self->{'pfilt'} $self->{'fasta_path'} > $self->{'infile'}";
    $self->{'cmd'}       = "$self->{'exe'} $self->{'infile'} > $self->{'outfile'}";
    
    return $self;
}

sub fasta_path
{
    print STDERR "**Seq2 fasta_path** \n";
    my ($self, $fasta_path) = @_;
    
    $self->{'fasta_path'} = $fasta_path if (defined($fasta_path));
    return $self->{'fasta_path'};
}

sub infile
{
    print STDERR "**seq2 infile** \n";
    my ($self, $infile) = @_;
    
    $self->{'infile'} = $infile if (defined($infile));
    return $self->{'infile'};
}

sub pfilt_cmd
{
    print STDERR "**Seq2 pfilt** \n";
    my ($self, $pfilt_cmd) = @_;
    
    $self->{'pfilt_cmd'} = $pfilt_cmd if (defined($pfilt_cmd));
    return $self->{'pfilt_cmd'};
}

sub run
{
    print STDERR "**Seq2 run** \n";
    my ($self, $error) = @_;
    $error = 0 unless ($error);
    my $infile = $self->infile();
    my $fasta_path = $self->fasta_path();
    
    if ($infile =~ /unmasked/)
    {
        copy($fasta_path, $infile);
    }
    else
    {
        print STDERR "Running pfilt for seq2mtx\n";
        $error += system($self->pfilt_cmd());
    }
    
    print STDERR "Running seq2mtx: ";
    $self->SUPER::run($error);
}



1;
