package AlignDB::Util;
use strict;
use warnings;
use autodie;

use Carp;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(firstidx all any uniq );
use Math::Combinatorics;
use Statistics::Descriptive;

use Bio::Seq;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::TCoffee;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::Align::DNAStatistics;

use AlignDB::IntSpan;

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA         = qw(Exporter);
%EXPORT_TAGS = (
    all => [
        qw{
            calc_gc_ratio pair_seq_stat single_indel_sites find_indel_set
            multi_align multi_align_matrix revcom seq_length
            average mean median variance stddev
            read_fasta write_fasta
            trim_pure_dash trim_head_tail trim_outgroup trim_complex_indel
            },
    ],
);
@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '1.0.1';

#----------------------------------------------------------#
# overide bioperl tempdir
#----------------------------------------------------------#
my $tmpdir;
{
    no warnings 'all';
    my $tempdir = sub {
        if ( defined $tmpdir ) {
            mkdir $tmpdir if !-e $tmpdir;
            return $tmpdir;
        }
        $tmpdir = 'tmp' . int( rand() * 10_000_000 );
        if ( $^O eq 'MSWin32' ) {
            $tmpdir = $ENV{TEMP} . '\\' . $tmpdir;
        }
        else {
            $tmpdir = '/tmp/' . $tmpdir;
        }
        mkdir $tmpdir if !-e $tmpdir;
        return $tmpdir;
    };
    *Bio::Tools::Run::Alignment::Clustalw::tempdir = $tempdir;
    *Bio::Tools::Run::Alignment::TCoffee::tempdir  = $tempdir;
    *Bio::Tools::Run::Alignment::Muscle::tempdir   = $tempdir;
    *Bio::Tools::Run::Alignment::MAFFT::tempdir    = $tempdir;
    use warnings;
}

#----------------------------------------------------------#
# every subroutines
#----------------------------------------------------------#

sub calc_gc_ratio {
    my @seqs = @_;

    for my $seq (@seqs) {
        _ref2str( \$seq );
    }

    my @ratios;
    for my $seq (@seqs) {

        # Count all four bases
        my $a_count = $seq =~ tr/Aa/Aa/;
        my $g_count = $seq =~ tr/Gg/Gg/;
        my $c_count = $seq =~ tr/Cc/Cc/;
        my $t_count = $seq =~ tr/Tt/Tt/;

        my $four_count = $a_count + $g_count + $c_count + $t_count;
        my $gc_count   = $g_count + $c_count;

        if ( $four_count == 0 ) {
            next;
        }
        else {
            my $gc_ratio = $gc_count / $four_count;
            push @ratios, $gc_ratio;
        }
    }

    return mean(@ratios);
}

sub single_indel_sites {
    my $seq = shift;

    _ref2str( \$seq );

    my $seq_legnth = length $seq;
    my @indel_sites;

    my $indel_offset = 0;
    my $indel_start  = 0;
    my $indel_end    = 0;
    for my $pos ( 1 .. $seq_legnth ) {
        my $base = substr( $seq, $pos - 1, 1 );
        if ( $base eq '-' ) {
            if ( $indel_offset == 0 ) {
                $indel_start = $pos;
            }
            $indel_offset++;
        }
        else {
            if ( $indel_offset != 0 ) {
                $indel_end = $pos - 1;
                my $indel_length = $indel_end - $indel_start + 1;

                push @indel_sites,
                    {
                    length => $indel_length,
                    start  => $indel_start,
                    end    => $indel_end,
                    };
            }
            $indel_offset = 0;
        }
    }
    if ( $indel_offset != 0 ) {
        $indel_end = $seq_legnth;
        my $indel_length = $indel_end - $indel_start + 1;

        push @indel_sites,
            {
            length => $indel_length,
            start  => $indel_start,
            end    => $indel_end,
            };
    }
    @indel_sites = sort { $a->{start} <=> $b->{start} } @indel_sites;

    return \@indel_sites;
}

sub find_indel_set {
    my $seq = shift;
    my $expand = shift || 0;

    my $seq_length = length $seq;

    my $indel_sites_ref = &single_indel_sites($seq);

    my $set = AlignDB::IntSpan->new();
    foreach (@$indel_sites_ref) {
        my $indel_runlist
            = ( $_->{start} - $expand ) . "-" . ( $_->{end} + $expand );
        $set = $set->add($indel_runlist);
    }

    $set = $set->intersect("1-$seq_length");

    return $set;
}

# clustalW align an array of sequences
sub multi_align {
    my ( $seqs_ref, $aln_prog ) = @_;

    my $seq_number = scalar @$seqs_ref;
    my ( @seqs, @seqs_obj );

    $aln_prog ||= 'clustalw';

    for ( my $i = 0; $i < $seq_number; $i++ ) {
        $seqs[$i] = $seqs_ref->[$i];
        $seqs[$i] =~ s/-//g;
        $seqs_obj[$i] = Bio::Seq->new(
            -display_id => "seq_$i",
            -seq        => $seqs[$i],
        );
    }

    my $aln_factory;
    if ( $aln_prog =~ /clus/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new( -verbose => -1 );
    }
    elsif ( $aln_prog =~ /t\_?cof/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::TCoffee->new( -verbose => -1 );
    }
    elsif ( $aln_prog =~ /musc/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::Muscle->new( -verbose => -1 );
    }
    elsif ( $aln_prog =~ /maff/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::MAFFT->new( -verbose => -1 );
    }
    else {
        confess "Provide clustalw, tcoffee, muscle or mafft" . " as alignment program names";
    }
    unless ( $aln_factory->executable ) {
        confess "Could not find the executable for $aln_prog\n";
    }

    my $dna_aln = $aln_factory->align( \@seqs_obj );

    # Does a s/$arg1/$arg2/ on the sequences. Useful for gap characters
    $dna_aln->map_chars( '\.', '-' );

    foreach my $seq ( $dna_aln->each_seq() ) {
        my $seq_id = $seq->display_id;
        $seq_id =~ /seq_(\d+)/;
        my $seq_sn  = $1;
        my $seq_seq = $seq->seq;
        $seqs[$seq_sn] = $seq->seq;
    }

    return \@seqs;
}

# quick calc distance matrix
sub multi_align_matrix {
    my $seqs_ref = shift;

    my $seq_number = scalar @{$seqs_ref};
    my ( @seqs, @seqs_obj );

    for ( my $i = 0; $i < $seq_number; $i++ ) {
        $seqs[$i] = $seqs_ref->[$i];
        $seqs[$i] =~ s/-//g;
        $seqs_obj[$i] = Bio::Seq->new(
            -display_id => "seq_$i",
            -seq        => $seqs[$i],
        );
    }

    my $aln_factory = Bio::Tools::Run::Alignment::MAFFT->new( -verbose => -1 );
    unless ( $aln_factory->executable ) {
        confess "Could not find the executable for mafft\n";
    }

    my $dna_aln = $aln_factory->align( \@seqs_obj );

    my $stats  = Bio::Align::DNAStatistics->new();
    my $matrix = $stats->distance(
        -align  => $dna_aln,
        -method => 'Jukes-Cantor'
    );

    return $matrix;
}

sub seq_length {
    my $seq = shift;

    _ref2str( \$seq );
    my $gaps = $seq =~ tr/-/-/;

    return length($seq) - $gaps;
}

sub revcom {
    my $seq = shift;

    _ref2str( \$seq );
    $seq =~ tr/ACGTMRWSYKVHDBNacgtmrwsykvhdbn-/TGCAKYWSRMBDHVNtgcakywsrmbdhvn-/;
    my $seq_rc = reverse $seq;

    return $seq_rc;
}

sub average {
    my @num = @_;
    if ( @num == 0 ) {
        return 0;
    }
    else {
        return mean(@num);
    }
}

# in situ convert reference of string to string
# For the sake of efficiency, the return value should be discarded
sub _ref2str {
    my $ref = shift;

    if ( ref $ref eq "REF" ) {
        $$ref = $$$ref;    # this is very weird, but it works
    }

    unless ( ref $ref eq "SCALAR" ) {
        carp "Wrong parameter passed\n";
    }

    return $ref;
}

sub mean {
    @_ = grep { defined $_ } @_;
    return unless @_;
    return $_[0] unless @_ > 1;
    return sum(@_) / scalar(@_);
}

sub median {
    @_ = grep { defined $_ } @_;
    return unless @_;
    return $_[0] unless @_ > 1;
    @_ = sort { $a <=> $b } @_;
    return $_[ $#_ / 2 ] if @_ & 1;
    my $mid = @_ / 2;
    return ( $_[ $mid - 1 ] + $_[$mid] ) / 2;
}

sub variance {
    @_ = grep { defined $_ } @_;
    return   unless @_;
    return 0 unless @_ > 1;
    my $mean = mean(@_);
    return sum( map { ( $_ - $mean )**2 } @_ ) / $#_;
}

sub stddev {
    @_ = grep { defined $_ } @_;
    return   unless @_;
    return 0 unless @_ > 1;
    return sqrt variance(@_);
}

sub read_fasta {
    my $filename = shift;

    my @lines = path($filename)->lines;

    my @seq_names;
    my %seqs;
    for my $line (@lines) {
        if ( $line =~ /^\>[\w:-]+/ ) {
            $line =~ s/\>//;
            chomp $line;
            push @seq_names, $line;
            $seqs{$line} = '';
        }
        elsif ( $line =~ /^[\w-]+/ ) {
            $line =~ s/[^\w-]//g;
            chomp $line;
            my $seq_name = $seq_names[-1];
            $seqs{$seq_name} .= $line;
        }
        else {    # Blank line, do nothing
        }
    }

    return ( \%seqs, \@seq_names );
}

sub write_fasta {
    my $filename   = shift;
    my $seq_of     = shift;
    my $seq_names  = shift;
    my $real_names = shift;

    open my $fh, ">", $filename;
    for my $i ( 0 .. @{$seq_names} - 1 ) {
        my $seq = $seq_of->{ $seq_names->[$i] };
        my $header;
        if ($real_names) {
            $header = $real_names->[$i];
        }
        else {
            $header = $seq_names->[$i];
        }

        print {$fh} ">" . $header . "\n";
        print {$fh} $seq . "\n";
    }
    close $fh;

    return;
}

#----------------------------#
# trim header and footer pure dash regions
#----------------------------#
sub trim_pure_dash {
    my $seq_of    = shift;
    my $seq_names = shift;

    # header indels
    while (1) {
        my @first_column;
        for ( @{$seq_names} ) {
            my $first_base = substr( $seq_of->{$_}, 0, 1 );
            push @first_column, $first_base;
        }
        if ( all { $_ eq '-' } @first_column ) {
            for ( @{$seq_names} ) {
                my $base = substr( $seq_of->{$_}, 0, 1, '' );
            }
        }
        else {
            last;
        }
    }

    # footer indels
    while (1) {
        my (@last_column);
        for ( @{$seq_names} ) {
            my $last_base = substr( $seq_of->{$_}, -1, 1 );
            push @last_column, $last_base;
        }
        if ( all { $_ eq '-' } @last_column ) {
            for ( @{$seq_names} ) {
                my $base = substr( $seq_of->{$_}, -1, 1, '' );
            }
        }
        else {
            last;
        }
    }

    return;
}

#----------------------------#
# trim head and tail indels
#----------------------------#
#  If head length set to 1, the first indel will be trimmed
#  Length set to 5 and the second indel will also be trimmed
#   GAAA--C
#   --AAAGC
#   GAAAAGC
sub trim_head_tail {
    my $seq_of      = shift;
    my $seq_names   = shift;
    my $head_length = shift;    # indels in this region will also be trimmed
    my $tail_length = shift;    # indels in this region will also be trimmed

    # default value means only trim indels starting at the first base
    $head_length = defined $head_length ? $head_length : 1;
    $tail_length = defined $tail_length ? $tail_length : 1;

    my $seq_number   = scalar @{$seq_names};
    my $align_length = length $seq_of->{ $seq_names->[0] };

    my $align_set = AlignDB::IntSpan->new("1-$align_length");
    my $indel_set = AlignDB::IntSpan->new;

    for my $n ( @{$seq_names} ) {
        my $seq_indel_set = find_indel_set( $seq_of->{$n} );
        $indel_set->merge($seq_indel_set);
    }

    # record bp chopped
    my %head_chopped = map { $_ => 0 } @{$seq_names};
    my %tail_chopped = map { $_ => 0 } @{$seq_names};

    # There're no indels at all
    return ( \%head_chopped, \%tail_chopped ) if $indel_set->is_empty;

    # head indel(s) to be trimmed
    my $head_set = AlignDB::IntSpan->new;
    $head_set->add_range( 1, $head_length );
    my $head_indel_set = $indel_set->find_islands($head_set);

    # head indels
    if ( $head_indel_set->is_not_empty ) {
        for my $i ( 1 .. $head_indel_set->max ) {
            my @column;
            for my $n ( @{$seq_names} ) {
                my $base = substr( $seq_of->{$n}, 0, 1, '' );
                if ( $base ne '-' ) {
                    $head_chopped{$n}++;
                }
            }
        }
    }

    # tail indel(s) to be trimmed
    my $tail_set = AlignDB::IntSpan->new;
    $tail_set->add_range( $align_length - $tail_length + 1, $align_length );
    my $tail_indel_set = $indel_set->find_islands($tail_set);

    # tail indels
    if ( $tail_indel_set->is_not_empty ) {
        for my $i ( $tail_indel_set->min .. $align_length ) {
            my @column;
            for my $n ( @{$seq_names} ) {
                my $base = substr( $seq_of->{$n}, -1, 1, '' );
                if ( $base ne '-' ) {
                    $tail_chopped{$n}++;
                }
            }
        }
    }

    return ( \%head_chopped, \%tail_chopped );
}

#----------------------------#
# trim outgroup only sequence
#----------------------------#
# if intersect is superset of union
#   target G----C
#   query  G----C
#   ref    GAAAAC
sub trim_outgroup {
    my $seq_of    = shift;
    my $seq_names = shift;

    # don't expand indel set here
    # last is outgroup
    my %indel_sets;
    for ( 0 .. @{$seq_names} - 2 ) {
        my $name = $seq_names->[$_];
        $indel_sets{$name} = find_indel_set( $seq_of->{$name} );
    }

    # find trim_region
    my $trim_region = AlignDB::IntSpan->new;

    my $union_set     = AlignDB::IntSpan::union( values %indel_sets );
    my $intersect_set = AlignDB::IntSpan::intersect( values %indel_sets );

    for my $span ( $union_set->runlists ) {
        if ( $intersect_set->superset($span) ) {
            $trim_region->add($span);
        }
    }

    # trim all segments in trim_region
    print " " x 4, "Delete trim region " . $trim_region->runlist . "\n"
        if $trim_region->is_not_empty;
    for my $span ( reverse $trim_region->spans ) {
        my $seg_start = $span->[0];
        my $seg_end   = $span->[1];
        for my $name ( @{$seq_names} ) {
            substr( $seq_of->{$name}, $seg_start - 1, $seg_end - $seg_start + 1, '' );
        }
    }

    return;
}

#----------------------------#
# record complex indels and ingroup indels
#----------------------------#
# if intersect is subset of union
#   tar GGA--C
#   que G----C
#   ref GGAGAC
sub trim_complex_indel {
    my $seq_of    = shift;
    my $seq_names = shift;

    my $ingroup_names = [ @{$seq_names} ];
    my $outgroup_name = pop @{$ingroup_names};

    my $complex_region = AlignDB::IntSpan->new;

    # don't expand indel set
    my %indel_sets;
    for ( @{$seq_names} ) {
        $indel_sets{$_} = find_indel_set( $seq_of->{$_} );
    }
    my $outgroup_indel_set = $indel_sets{$outgroup_name};
    delete $indel_sets{$outgroup_name};

    # all ingroup intersect sets are complex region after remove uniform ingroup
    #   indels
    my $union_set     = AlignDB::IntSpan::union( values %indel_sets );
    my $intersect_set = AlignDB::IntSpan::intersect( values %indel_sets );

    print " " x 4, "Delete complex trim region " . $intersect_set->runlist . "\n"
        if $intersect_set->is_not_empty;
    for ( reverse $intersect_set->spans ) {
        my $seg_start = $_->[0];
        my $seg_end   = $_->[1];

        # trim sequence
        for ( @{$seq_names} ) {
            substr( $seq_of->{$_}, $seg_start - 1, $seg_end - $seg_start + 1, '' );
        }

        # add to complex_region
        for my $span ( $union_set->runlists ) {
            my $sub_union_set = AlignDB::IntSpan->new($span);
            if ( $sub_union_set->superset("$seg_start-$seg_end") ) {
                $complex_region->merge($sub_union_set);
            }
        }

        # modify all related set
        $union_set = $union_set->banish_span( $seg_start, $seg_end );
        for ( @{$ingroup_names} ) {
            $indel_sets{$_}
                = $indel_sets{$_}->banish_span( $seg_start, $seg_end );
        }
        $outgroup_indel_set->banish_span( $seg_start, $seg_end );
        $complex_region = $complex_region->banish_span( $seg_start, $seg_end );
    }

    # add ingroup-outgroup complex indels to complex_region
    for my $name ( @{$ingroup_names} ) {
        my $outgroup_intersect_set = $outgroup_indel_set->intersect( $indel_sets{$name} );
        for my $out_span ( $outgroup_intersect_set->runlists ) {
            for my $union_span ( $union_set->runlists ) {
                my $sub_union_set = AlignDB::IntSpan->new($union_span);

                # union_set > intersect_set
                if ( $sub_union_set->larger_than($out_span) ) {
                    $complex_region->merge($sub_union_set);
                }
            }
        }
    }

    return $complex_region->runlist;
}

1;

__END__


=pod

=encoding UTF-8

=head1 NAME

AlignDB::Util - Misc functions for AlignDB

=head1 SYNOPSIS

    use AlignDB::Util qw(:all);

=head1 AUTHOR

Qiang Wang <wang-q@outlook.com>

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2008- by Qiang Wang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
