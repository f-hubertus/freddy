#!/usr/bin/perl -w

package Tool;

use strict;
print "***Start Tool-Blast*** \n";


sub new
{
    my ($class) = @_;
    
    my $self = {
                 'cmd' => "",
                 'err' => 0
               };
    
    bless $self, $class;
    return $self;
}

sub cmd
{
    my ($self, $cmd) = @_;
    
    $self->{'cmd'} = $cmd if (defined($cmd));
    return $self->{'cmd'};
}

sub err
{
    my ($self, $err) = @_;
    
    $self->{'err'} = $err if (defined($err));
    return $self->{'err'};
}

sub return_type
{
    my ($self, $mask, $iter) = @_;
    
    return "" unless (exists($self->{'mask'}) && exists($self->{'iter'}));
    return $mask . $iter;
}

sub run
{
    my ($self, $error) = @_;
    $error = 0 unless ($error);
    
    print STDERR $self->cmd() . "**Tool-running** \n";
    $error += system($self->cmd());
    $self->err($error);
}



1;
