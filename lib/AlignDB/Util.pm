package AlignDB::Util;

# ABSTRACT: Misc functions for AlignDB

use strict;
use warnings;
use Carp;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Slurp;
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
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            calc_gc_ratio pair_seq_stat multi_seq_stat pair_snp_sites multi_snp_site
            single_indel_sites pair_indel_sites find_indel_set ref_indel_type ref_pair_D multi_align
            multi_align_matrix revcom seq_length average mean median variance stddev
            read_fasta write_fasta trim_pure_dash trim_head_tail trim_outgroup trim_complex_indel
            realign_quick decode_header encode_header
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

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

#----------------------------#
# Below are every approaches I have tried, just for memory
#----------------------------#

# AlignDB::Util::clustal_align use Bio::Tools::Run::Alignment::Clustalw, which
#   use File::Temp to generate temp files. AND File::Temp use File::Spec to
#   get temporary dir.
# So we change File::Spec::Unix::_tmpdir to get tmpdir relocated and now two
#   or more ref_outgroup.pl will not affect each other.
# PS: ingore the "Subroutine _tmpdir redefined" warning.
#package File::Spec::Unix;
#my $tmpdir;
#sub _tmpdir {
#    return $tmpdir if defined $tmpdir;
#    my $tmpdir = 'tmp' . int(rand() * 1000000000);
#    if ( $^O eq 'MSWin32' ) {
#        $tmpdir = 'c:\\' . $tmpdir ;
#    }
#    else {
#        $tmpdir = '/tmp/' . $tmpdir ;
#    }
#    mkdir $tmpdir if ! -e $tmpdir;
#    return $tmpdir;
#}

# Damn! Bio::Tools::Run:: also define a tempdir method. We should only
#   override it.
# We are really happy now!
# I failed to use Sub::Override, glob or even Symbol::Glob to replace this ugly
#   approach, is there a clean solution?
#package Bio::Tools::Run::WrapperBase;
#
#my $tmpdir;
#
#sub tempdir {
#    if ( defined $tmpdir ) {
#        mkdir $tmpdir if !-e $tmpdir;
#        return $tmpdir;
#    }
#    $tmpdir = 'tmp' . int( rand() * 1000000 );
#    if ( $^O eq 'MSWin32' ) {
#        $tmpdir = $ENV{TEMP} . '\\' . $tmpdir;
#    }
#    else {
#        $tmpdir = '/tmp/' . $tmpdir;
#    }
#    mkdir $tmpdir if !-e $tmpdir;
#    return $tmpdir;
#}

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

sub pair_seq_stat {
    my ( $first_seq, $second_seq ) = @_;

    _ref2str( \$first_seq );
    _ref2str( \$second_seq );

    my $seq_legnth = length $first_seq;

    # For every positions, search for polymorphism_site
    my ( $number_of_comparable_bases, $number_of_identities,
        $number_of_differences, $number_of_gaps, $number_of_n, $number_of_align_error, )
        = (0) x 6;
    for my $pos ( 1 .. $seq_legnth ) {
        my @nt_pair = ();
        foreach ( $first_seq, $second_seq ) {
            my $nt = substr( $_, $pos - 1, 1 );
            unless ( defined $nt ) {
                $nt = '?';
            }
            push @nt_pair, $nt;
        }
        if ( $nt_pair[0] =~ /[agct]/i ) {
            if ( $nt_pair[1] =~ /[agct]/i ) {
                if ( $nt_pair[0] ne $nt_pair[1] ) {
                    $number_of_differences++;
                }
                else {
                    $number_of_identities++;
                }
                $number_of_comparable_bases++;
            }
            elsif ( $nt_pair[1] =~ /\-/i ) {
                $number_of_gaps++;
            }
            else {
                $number_of_n++;
            }
        }
        elsif ( $nt_pair[0] =~ /\-/i ) {
            if ( $nt_pair[1] =~ /[agct]/i ) {
                $number_of_gaps++;
            }
            else {
                $number_of_align_error++;
            }
        }
        else {
            if ( $nt_pair[1] =~ /[agct]/i ) {
                $number_of_n++;
            }
            else {
                $number_of_align_error++;
            }
        }
    }
    if ( $number_of_comparable_bases == 0 ) {
        print Dump(
            {   first_seq  => $first_seq,
                second_seq => $second_seq,
            }
        );
        carp "number_of_comparable_bases == 0!!\n";
        return [
            $seq_legnth,            $number_of_comparable_bases, $number_of_identities,
            $number_of_differences, $number_of_gaps,             $number_of_n,
            $seq_legnth,            'NULL',                      'NULL',
            'NULL',
        ];
    }
    my $pi = $number_of_differences / $number_of_comparable_bases;

    my $first_seq_gc = calc_gc_ratio($first_seq);
    my $average_gc = calc_gc_ratio( $first_seq, $second_seq );

    return [
        $seq_legnth,            $number_of_comparable_bases, $number_of_identities,
        $number_of_differences, $number_of_gaps,             $number_of_n,
        $number_of_align_error, $pi,                         $first_seq_gc,
        $average_gc,
    ];
}

sub multi_seq_stat {
    my (@seqs) = @_;

    for my $seq (@seqs) {
        _ref2str( \$seq );
    }

    my $seq_legnth = length $seqs[0];

    # For every positions, search for polymorphism_site
    my ( $number_of_comparable_bases, $number_of_identities,
        $number_of_differences, $number_of_gaps, $number_of_n, $number_of_align_error, )
        = (0) x 6;
    for my $pos ( 1 .. $seq_legnth ) {
        my @bases = ();
        for (@seqs) {
            my $nt = substr( $_, $pos - 1, 1 );
            push @bases, $nt;
        }
        @bases = uniq(@bases);

        if ( all { $_ =~ /[agct]/i } @bases ) {
            $number_of_comparable_bases++;
            if ( all { $_ eq $bases[0] } @bases ) {
                $number_of_identities++;
            }
            else {
                $number_of_differences++;
            }
        }
        elsif ( any { $_ eq '-' } @bases ) {
            $number_of_gaps++;
        }
        else {
            $number_of_n++;
        }
    }
    if ( $number_of_comparable_bases == 0 ) {
        print Dump { seqs => \@seqs, };
        carp "number_of_comparable_bases == 0!!\n";
        return [
            $seq_legnth, $number_of_comparable_bases, $number_of_identities, $number_of_differences,
            $number_of_gaps, $number_of_n, $seq_legnth, undef, undef, undef,
        ];
    }

    my $combinat = Math::Combinatorics->new(
        count => 2,
        data  => [ 0 .. @seqs - 1 ],
    );
    my @all_pi;
    while ( my ( $idx1, $idx2 ) = $combinat->next_combination ) {
        my $stats = pair_seq_stat( $seqs[$idx1], $seqs[$idx2] );
        my $pi = $stats->[7];
        push @all_pi, $pi;
    }
    my $pi = mean(@all_pi);

    my $target_gc  = calc_gc_ratio( $seqs[0] );
    my $average_gc = calc_gc_ratio(@seqs);

    return [
        $seq_legnth,            $number_of_comparable_bases, $number_of_identities,
        $number_of_differences, $number_of_gaps,             $number_of_n,
        $number_of_align_error, $pi,                         $target_gc,
        $average_gc,
    ];
}

sub pair_snp_sites {
    my ( $first_seq, $second_seq ) = @_;

    _ref2str( \$first_seq );
    _ref2str( \$second_seq );

    my $seq_legnth = length $first_seq;
    my %snp_sites;

    for my $pos ( 1 .. $seq_legnth ) {
        my @nt_pair = ();
        foreach ( $first_seq, $second_seq ) {
            my $nt = substr $_, $pos - 1, 1;
            unless ( defined $nt ) {
                $nt = '?';
            }
            push @nt_pair, $nt;
        }
        if ( $nt_pair[0] =~ /[agct]/i ) {
            if ( $nt_pair[1] =~ /[agct]/i ) {
                if ( $nt_pair[0] ne $nt_pair[1] ) {
                    $snp_sites{$pos}->{target_base} = $nt_pair[0];
                    $snp_sites{$pos}->{query_base}  = $nt_pair[1];
                }
            }
        }
    }

    return \%snp_sites;
}

sub multi_snp_site {
    my (@seqs) = @_;

    for my $seq (@seqs) {
        _ref2str( \$seq );
    }

    my $seq_legnth = length $seqs[0];
    my %snp_sites;
    for my $pos ( 1 .. $seq_legnth ) {
        my @bases = ();
        foreach (@seqs) {
            my $nt = substr( $_, $pos - 1, 1 );
            push @bases, $nt;
        }
        if ( all { $_ =~ /[agct]/i } @bases ) {
            if ( any { $_ ne $bases[0] } @bases ) {
                $snp_sites{$pos} = \@bases;
            }
        }
        else {
            next;
        }
    }

    return \%snp_sites;
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

sub pair_indel_sites {
    my ( $first_seq, $second_seq, $defined_indel ) = @_;

    _ref2str( \$first_seq );
    _ref2str( \$second_seq );

    my $indel_set = AlignDB::IntSpan->new();
    if ( defined $defined_indel ) {
        $indel_set->merge($defined_indel);
    }
    else {
        $indel_set->merge( &find_indel_set($first_seq) );
        $indel_set->merge( &find_indel_set($second_seq) );
    }

    my @indel_sites;
    foreach my $span ( $indel_set->spans() ) {
        my $indel_start  = $span->[0];
        my $indel_end    = $span->[1];
        my $indel_length = $indel_end - $indel_start + 1;

        my $first_indel_seq  = substr( $first_seq,  $indel_start - 1, $indel_length );
        my $second_indel_seq = substr( $second_seq, $indel_start - 1, $indel_length );

        # $indel_insert:
        #   'N': means indel occured in other place
        #   'D': means deletion relative to first seq
        #   'I': means insertion relative to first seq
        # $indel_seq:
        #   'N': use $first_indel_seq as $indel_seq
        #   'D': use $first_indel_seq as $indel_seq
        #   'I': use $second_indel_seq as $indel_seq
        my $indel_insert;
        my $indel_seq;
        if ( $first_seq eq $second_seq ) {
            $indel_insert = 'N';
            $indel_seq    = $first_indel_seq;
        }
        elsif ( $first_seq !~ /\-/ and $second_seq !~ /\-/ ) {
            $indel_insert = 'N';
            $indel_seq    = $first_indel_seq;
        }
        else {
            my $first_gap  = $first_indel_seq =~ tr/-/-/;
            my $second_gap = $second_indel_seq =~ tr/-/-/;
            if ( $first_gap < $second_gap ) {
                $indel_insert = 'D';
                $indel_seq    = $first_indel_seq;
            }
            elsif ( $first_gap > $second_gap ) {
                $indel_insert = 'I';
                $indel_seq    = $second_indel_seq;
            }
            else {
                $indel_insert = 'N';
                $indel_seq    = $first_indel_seq;
            }
        }

        my $indel_gc = &calc_gc_ratio($indel_seq);
        push @indel_sites,
            {
            insert => $indel_insert,
            length => $indel_length,
            start  => $indel_start,
            end    => $indel_end,
            seq    => $indel_seq,
            gc     => $indel_gc,
            };
    }

    my $seq_legnth = length $first_seq;

    my $anterior_indel_end = 0;
    for my $i ( 0 .. scalar @indel_sites - 1 ) {
        my $current_indel_start = $indel_sites[$i]->{start};
        my $current_left_extand = $current_indel_start - 1 - $anterior_indel_end;
        my $current_indel_end   = $indel_sites[$i]->{end};
        $anterior_indel_end = $current_indel_end;
        $indel_sites[$i]->{left_extand} = $current_left_extand;
    }

    my $posterior_indel_start = $seq_legnth + 1;
    for my $i ( reverse( 0 .. scalar @indel_sites - 1 ) ) {
        my $current_indel_end    = $indel_sites[$i]->{end};
        my $current_right_extand = $posterior_indel_start - $current_indel_end - 1;
        my $current_indel_start  = $indel_sites[$i]->{start};
        $posterior_indel_start = $current_indel_start;
        $indel_sites[$i]->{right_extand} = $current_right_extand;
    }

    return \@indel_sites;
}

sub ref_indel_type {
    my $ref_str    = shift;
    my $target_str = shift;
    my $query_str  = shift;

    my $rindel_set = &find_indel_set($ref_str);
    my $tindel_set = &find_indel_set($target_str);
    my $qindel_set = &find_indel_set($query_str);

    #print Dump {
    #    "1rindel" => "$ref_str $rindel_set",
    #    "2tindel" => "$target_str $tindel_set",
    #    "3qindel" => "$query_str $qindel_set",
    #};

    # $indel_occured:
    #   'C': complex indel
    #   'T': occured in target seq
    #   'Q': occured in query seq
    #   'N': occured in other place, noindel
    # $indel_type:
    #   'C': complex indel
    #   'D': deletion
    #   'I': insertion
    #   'N': noindel
    my $indel_occured = "NULL";
    my $indel_type    = "NULL";

    my $in_intersect = $tindel_set->intersect($qindel_set);
    my $in_union     = $tindel_set->union($qindel_set);

    if ( $in_intersect->equal($in_union) ) {
        $indel_occured = "N";
        $indel_type    = "N";
    }
    elsif ( $rindel_set->equal($in_union) ) {
        if ( $tindel_set->larger_than($qindel_set) ) {
            $indel_occured = "Q";
            $indel_type    = "I";
        }
        elsif ( $tindel_set->smaller_than($qindel_set) ) {
            $indel_occured = "T";
            $indel_type    = "I";
        }
        else {
            $indel_occured = "C";
            $indel_type    = "C";
        }
    }
    elsif ( $rindel_set->equal($in_intersect) ) {
        if ( $tindel_set->larger_than($qindel_set) ) {
            $indel_occured = "T";
            $indel_type    = "D";
        }
        elsif ( $tindel_set->smaller_than($qindel_set) ) {
            $indel_occured = "Q";
            $indel_type    = "D";
        }
        else {
            $indel_occured = "C";
            $indel_type    = "C";
        }
    }
    elsif ( $rindel_set->larger_than($in_union) ) {
        $indel_occured = "C";
        $indel_type    = "C";
    }
    elsif ( $rindel_set->smaller_than($in_union) ) {
        $indel_occured = "C";
        $indel_type    = "C";
    }
    else {
        $indel_occured = "C";
        $indel_type    = "C";
    }

    return ( $indel_occured, $indel_type );
}

##################################################
# Usage      : my ( $d1, $d2, $dc ) = &ref_pair_D(
#            :     $ref_seq,
#            :     $first_seq,
#            :     $second_seq
#            : );
# Purpose    : Split D value to D1 (substitutions in first_seq),
#            : D2( substitutions in second_seq) and
#            : Dcomplex (substitutions can't be referred)
# Returns    : ( $d1, $d2, $dc )
# Parameters : $ref_seq, $first_seq, $second_seq
# Throws     : no exceptions
# Comments   : If use $indel_side_seq and $second_seq as $first_seq
#            : and $second_seq, $d1 and $d2 will be $di and $dn
#            : Modified from Dr. Zhu's codes
# See Also   : none
sub ref_pair_D {
    my $ref_str    = shift;
    my $first_str  = shift;
    my $second_str = shift;

    my ( $d1, $d2, $dc ) = (0) x 3;
    my $length = length $ref_str;

    return ( $d1, $d2, $dc ) if $length == 0;

    for my $i ( 0 .. $length - 1 ) {
        my $first_base  = substr $first_str,  $i, 1;
        my $second_base = substr $second_str, $i, 1;
        my $ref_base    = substr $ref_str,    $i, 1;
        if ( $first_base ne $second_base ) {
            if (   $first_base =~ /[ATCG]/i
                && $second_base =~ /[ATCG]/i
                && $ref_base =~ /[ATCG]/i )
            {
                if ( $second_base eq $ref_base ) {
                    $d1++;
                }
                elsif ( $first_base eq $ref_base ) {
                    $d2++;
                }
                else {
                    $dc++;
                }
            }
            else {
                $dc++;
            }
        }
    }

    for ( $d1, $d2, $dc ) {
        $_ /= $length;
    }

    return ( $d1, $d2, $dc );
}

# clustalW align an array of sequences
sub clustal_align {
    my ( $seqs_ref, $more_indel ) = @_;

    my $seq_number = scalar @$seqs_ref;
    my ( @seqs, @seqs_obj );

    for ( my $i = 0; $i < $seq_number; $i++ ) {
        $seqs[$i] = $seqs_ref->[$i];
        $seqs[$i] =~ s/-//g;
        $seqs_obj[$i] = Bio::Seq->new(
            -display_id => "seq_$i",
            -seq        => $seqs[$i],
        );
    }

    # Build a clustalw alignment factory
    # Key parameters are gap-related ones
    #   pairwise alignments
    #     -pairgap=n   Gap penalty.
    #     -pwgapopen=f Gap opening penalty.     default: 15.00
    #     -pwgapext=f  Gap extension penalty.   default: 6.66
    #   multiple alignments
    #     -gapdist=n   Gap separation penalty range.
    #     -gapopen=f   Gap opening penalty.     default: 15.0
    #     -gapext=f    Gap extension penalty.   default: 6.66
    my @params;

    if ($more_indel) {
        @params = (
            'pwgapopen' => 5,
            'pwgapext'  => 2.22,
            'gapopen'   => 5,
            'gapext'    => 2.22,
        );
    }
    my $aln_factory;
    $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new( @params, -verbose => -1 );
    unless ( $aln_factory->executable ) {
        die "Could not find the executable for ClustalW\n";
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

sub factorial_of {
    my $n = shift;
    return unless $n >= 0 and $n == int($n);

    my $f;

    for ( $f = 1; $n > 0; $n-- ) {
        $f *= $n;
    }

    return $f;
}

#sub combi_k_n {
#    my ( $k, $n ) = @_;
#    my $com = factorial_of($n) / factorial_of($k) / factorial_of( $n - $k );
#    return $com;
#}

## Usage
##
##my @comb = (-1);
##
##while (&enumComb(10, 3, \@comb)) {
##    print "@comb\n";
##}
#sub enumComb {
#    my ( $n, $k, $comb ) = @_;
#    my $i;
#    if ( $comb->[0] < 0 ) {
#        for ( $i = 0; $i < $k; $i++ ) {
#            $comb->[$i] = $i;
#        }
#        return 1;
#    }
#    else {
#        for ( $i = $k - 1; $i >= 0 && $comb->[$i] >= $n - $k + $i; $i-- ) { }
#        if ( $i >= 0 ) {
#            $comb->[$i]++;
#            my $m;
#            for ( $m = $i + 1; $m < $k; $m++ ) {
#                $comb->[$m] = $comb->[ $m - 1 ] + 1;
#            }
#            return 1;
#        }
#        else {
#            return 0;
#        }
#    }
#}

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

    my @lines = read_file($filename);

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

#----------------------------#
# realign indel_flank region
#----------------------------#
sub realign_quick {
    my $seq_of    = shift;
    my $seq_names = shift;
    my $opt       = shift;

    if ( !exists $opt->{indel_expand} ) {
        $opt->{indel_expand} = 50;
    }
    if ( !exists $opt->{indel_join} ) {
        $opt->{indel_join} = 50;
    }
    if ( !exists $opt->{aln_prog} ) {
        $opt->{aln_prog} = 'clustalw';
    }

    # use AlignDB::IntSpan to find nearby indels
    #   expand indel by a range of $indel_expand
    my %indel_sets;
    for (@$seq_names) {
        $indel_sets{$_} = find_indel_set( $seq_of->{$_}, $opt->{indel_expand} );
    }

    my $realign_region = AlignDB::IntSpan->new;
    my $combinat       = Math::Combinatorics->new(
        count => 2,
        data  => $seq_names,
    );
    while ( my @combo = $combinat->next_combination ) {
        my $intersect_set = AlignDB::IntSpan->new;
        my $union_set     = AlignDB::IntSpan->new;
        $intersect_set
            = $indel_sets{ $combo[0] }->intersect( $indel_sets{ $combo[1] } );
        $union_set
            = $indel_sets{ $combo[0] }->union( $indel_sets{ $combo[1] } );

        for my $span ( $union_set->runlists ) {
            my $flag_set = $intersect_set->intersect($span);
            if ( $flag_set->is_not_empty ) {
                $realign_region->add($span);
            }
        }
    }

    # join adjacent realign regions
    $realign_region = $realign_region->join_span( $opt->{indel_join} );

    # realign all segments in realign_region
    my @realign_region_spans = $realign_region->spans;
    for ( reverse @realign_region_spans ) {
        my $seg_start = $_->[0];
        my $seg_end   = $_->[1];
        my @segments;
        for (@$seq_names) {
            my $seg = substr( $seq_of->{$_}, $seg_start - 1, $seg_end - $seg_start + 1 );
            push @segments, $seg;
        }

        my $realigned_segments = multi_align( \@segments, $opt->{aln_prog} );

        for (@$seq_names) {
            my $seg = shift @{$realigned_segments};
            $seg = uc $seg;
            substr( $seq_of->{$_}, $seg_start - 1, $seg_end - $seg_start + 1, $seg );
        }
    }

    return;
}

sub decode_header {
    my $header = shift;

    # S288C.chrI(+):27070-29557|species=S288C
    my $head_qr = qr{
                ([\w_]+)?           # name
                [\.]?               # spacer
                ((?:chr)?[\w-]+)    # chr name
                (?:\((.+)\))?       # strand
                [\:]                # spacer
                (\d+)               # chr start
                [\_\-]              # spacer
                (\d+)               # chr end
            }xi;

    my $info = {};

    $header =~ $head_qr;
    my $name     = $1;
    my $chr_name = $2;

    if ( defined $name ) {
        $info = {
            chr_name   => $2,
            chr_strand => $3,
            chr_start  => $4,
            chr_end    => $5,
        };
        if ( $info->{chr_strand} eq '1' ) {
            $info->{chr_strand} = '+';
        }
        elsif ( $info->{chr_strand} eq '-1' ) {
            $info->{chr_strand} = '-';
        }
    }
    elsif ( defined $chr_name ) {
        $name = $header;
        $info = {
            chr_name   => $2,
            chr_strand => $3,
            chr_start  => $4,
            chr_end    => $5,
        };
        if ( !defined $info->{chr_strand} ) {
            $info->{chr_strand} = '+';
        }
        elsif ( $info->{chr_strand} eq '1' ) {
            $info->{chr_strand} = '+';
        }
        elsif ( $info->{chr_strand} eq '-1' ) {
            $info->{chr_strand} = '-';
        }
    }
    else {
        $name = $header;
        $info = {
            chr_name   => 'chrUn',
            chr_strand => '+',
            chr_start  => undef,
            chr_end    => undef,
        };
    }
    $info->{name} = $name;

    # additional keys
    if ( $header =~ /\|(.+)/ ) {
        my @parts = grep {defined} split /;/, $1;
        for my $part (@parts) {
            my ( $key, $value ) = split /=/, $part;
            if ( defined $key and defined $value ) {
                $info->{$key} = $value;
            }
        }
    }

    return $info;
}

sub encode_header {
    my $info = shift;

    my $header;
    $header .= $info->{name};
    $header .= "." . $info->{chr_name};
    $header .= "(" . $info->{chr_strand} . ")";
    $header .= ":" . $info->{chr_start};
    $header .= "-" . $info->{chr_end};

    # additional keys
    my %essential = map { $_ => 1 } qw{name chr_name chr_strand chr_start chr_end seq full_seq};
    my @parts;
    for my $key ( sort keys %{$info} ) {
        if ( !$essential{$key} ) {
            push @parts, $key . "=" . $info->{$key};
        }
    }
    if (@parts) {
        my $additional = join ";", @parts;
        $header .= "|" . $additional;
    }

    return $header;
}

1;

__END__

=head1 SYNOPSIS

    use AlignDB::Util qw(:all);

