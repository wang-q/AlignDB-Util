package AlignDB::Util;

# ABSTRACT: Misc functions for AlignDB

use strict;
use warnings;
use Carp;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Slurp;
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

use AlignDB::IntSpan;

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            calc_gc_ratio
            pair_seq_stat         multi_seq_stat
            pair_snp_sites        multi_snp_site
            single_indel_sites    pair_indel_sites
            find_indel_set
            ref_indel_type        ref_pair_D
            clustal_align         multi_align
            random_sampling
            combi_k_n
            enumComb              k_nuc_permu
            k_nuc_count           k_nuc_incr
            revcom                seq_length
            average
            sampling_with_replacement
            random_number         stat_result
            mean                  median
            variance              stddev
            read_fasta
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
        $seq = uc $seq;
    }

    my @ratios;
    for my $seq (@seqs) {

        # Count all four bases
        my $a_count = $seq =~ tr/A/A/;
        my $g_count = $seq =~ tr/G/G/;
        my $c_count = $seq =~ tr/C/C/;
        my $t_count = $seq =~ tr/T/T/;

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
        $number_of_differences, $number_of_gaps, $number_of_n,
        $number_of_align_error, )
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
            $seq_legnth,           $number_of_comparable_bases,
            $number_of_identities, $number_of_differences,
            $number_of_gaps,       $number_of_n,
            $seq_legnth,           'NULL',
            'NULL',                'NULL',
        ];
    }
    my $pi = $number_of_differences / $number_of_comparable_bases;

    my $first_seq_gc = calc_gc_ratio($first_seq);
    my $average_gc = calc_gc_ratio( $first_seq, $second_seq );

    return [
        $seq_legnth,            $number_of_comparable_bases,
        $number_of_identities,  $number_of_differences,
        $number_of_gaps,        $number_of_n,
        $number_of_align_error, $pi,
        $first_seq_gc,          $average_gc,
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
        $number_of_differences, $number_of_gaps, $number_of_n,
        $number_of_align_error, )
        = (0) x 6;
    for my $pos ( 1 .. $seq_legnth ) {
        my @bases = ();
        foreach (@seqs) {
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
            $seq_legnth,           $number_of_comparable_bases,
            $number_of_identities, $number_of_differences,
            $number_of_gaps,       $number_of_n,
            $seq_legnth,           'NULL',
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
        $seq_legnth,            $number_of_comparable_bases,
        $number_of_identities,  $number_of_differences,
        $number_of_gaps,        $number_of_n,
        $number_of_align_error, $pi,
        $target_gc,             $average_gc,
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

        my $first_indel_seq
            = substr( $first_seq, $indel_start - 1, $indel_length );
        my $second_indel_seq
            = substr( $second_seq, $indel_start - 1, $indel_length );

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
            my $first_gap  = $first_indel_seq  =~ tr/-/-/;
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
        my $current_left_extand
            = $current_indel_start - 1 - $anterior_indel_end;
        my $current_indel_end = $indel_sites[$i]->{end};
        $anterior_indel_end = $current_indel_end;
        $indel_sites[$i]->{left_extand} = $current_left_extand;
    }

    my $posterior_indel_start = $seq_legnth + 1;
    for my $i ( reverse( 0 .. scalar @indel_sites - 1 ) ) {
        my $current_indel_end = $indel_sites[$i]->{end};
        my $current_right_extand
            = $posterior_indel_start - $current_indel_end - 1;
        my $current_indel_start = $indel_sites[$i]->{start};
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
    $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
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
        $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new();
    }
    elsif ( $aln_prog =~ /t\_?cof/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::TCoffee->new();
    }
    elsif ( $aln_prog =~ /musc/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::Muscle->new();
    }
    elsif ( $aln_prog =~ /maff/i ) {
        $aln_factory = Bio::Tools::Run::Alignment::MAFFT->new();
    }
    else {
        confess "Provide clustalw, tcoffee, muscle or mafft"
            . " as alignment program names";
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

# To select n records at random from a set of N, where 0 < n <= N
# return an array containing 0 .. N - 1
# Algorithm S (Selection sampling technique)
# TAOCP Vol2 3.4.2
sub random_sampling {
    my ( $N, $n ) = @_;

    my $t = 0;    # t is the total number of input records we have dealt with
    my $m = 0;    # m represents the number of records selected so far

    my @samples;

    while (1) {
        my $U = rand();
        if ( ( $N - $t ) * $U >= $n - $m ) {
            $t++;
        }
        else {
            push @samples, $t;
            $m++;
            $t++;
            last if ( $m >= $n );
        }
    }
    return @samples;
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

sub combi_k_n {
    my ( $k, $n ) = @_;
    my $com = factorial_of($n) / factorial_of($k) / factorial_of( $n - $k );
    return $com;
}

# Usage
#
#my @comb = (-1);
#
#while (&enumComb(10, 3, \@comb)) {
#    print "@comb\n";
#}
sub enumComb {
    my ( $n, $k, $comb ) = @_;
    my $i;
    if ( $comb->[0] < 0 ) {
        for ( $i = 0; $i < $k; $i++ ) {
            $comb->[$i] = $i;
        }
        return 1;
    }
    else {
        for ( $i = $k - 1; $i >= 0 && $comb->[$i] >= $n - $k + $i; $i-- ) { }
        if ( $i >= 0 ) {
            $comb->[$i]++;
            my $m;
            for ( $m = $i + 1; $m < $k; $m++ ) {
                $comb->[$m] = $comb->[ $m - 1 ] + 1;
            }
            return 1;
        }
        else {
            return 0;
        }
    }
}

sub k_nuc_permu {
    my $k         = shift;
    my @alphabets = qw{A C G T};

    my %table;
    $table{$_} = '' foreach @alphabets;
    foreach ( 2 .. $k ) {
        foreach my $current_key ( keys %table ) {
            $table{ $current_key . $_ } = '' foreach @alphabets;
            delete $table{$current_key};
        }
    }

    return sort keys %table;
}

sub k_nuc_count {
    my $seq_ref = shift;
    my $k       = shift;

    my $seq_length = length $$seq_ref;
    my %table;

    foreach ( 0 .. $seq_length - $k ) {
        $table{ substr( $$seq_ref, $_, $k ) }++;
    }

    return %table;
}

sub k_nuc_incr {
    my $seq_ref  = shift;
    my $k        = shift;
    my $hash_ref = shift;

    my $seq_length = length $$seq_ref;

    foreach ( 0 .. $seq_length - $k ) {
        $hash_ref->{ substr( $$seq_ref, $_, $k ) }++;
    }
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
    $seq
        =~ tr/ACGTMRWSYKVHDBNacgtmrwsykvhdbn-/TGCAKYSWRMBDHVNtgcakyswrmbdhvn-/;
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

sub sampling_with_replacement {
    my $data = shift;

    my @sample;
    my $size = scalar @$data;
    for ( 1 .. $size ) {
        my $random_index = random_number($size);
        push @sample, $data->[$random_index];
    }

    return \@sample;
}

sub random_number {
    my $max_int = shift;
    my $random  = int( rand($max_int) );
    return $random;
}

sub stat_result {
    my $data = shift;

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@$data);

    my $count    = $stat->count();
    my $mean     = $stat->mean();
    my $variance = $stat->variance();
    my $stddev   = $stat->standard_deviation();

    # standard_error =  standard_deviation / sqrt(sample_size)
    my $stderr = $stddev / sqrt($count);

    my $median        = $stat->median();
    my $percentile_5  = $stat->percentile(5);
    my $percentile_95 = $stat->percentile(95);

    my $result = {
        a_mean          => $mean,
        b_stderr        => $stderr,
        c_median        => $median,
        d_percentile_5  => $percentile_5,
        e_percentile_95 => $percentile_95,
    };

    return $result;
}

sub mean {
    return unless @_;
    return $_[0] unless @_ > 1;
    return sum(@_) / scalar(@_);
}

sub median {
    return unless @_;
    return $_[0] unless @_ > 1;
    @_ = sort { $a <=> $b } @_;
    return $_[ $#_ / 2 ] if @_ & 1;
    my $mid = @_ / 2;
    return ( $_[ $mid - 1 ] + $_[$mid] ) / 2;
}

sub variance {
    return   unless @_;
    return 0 unless @_ > 1;
    my $mean = mean(@_);
    return sum( map { ( $_ - $mean )**2 } @_ ) / $#_;
}

sub stddev {
    return   unless @_;
    return 0 unless @_ > 1;
    return sqrt variance(@_);
}

sub read_fasta {
    my $filename = shift;

    my @lines = read_file($filename);

    my @seq_names;
    my %seqs;
    foreach my $line (@lines) {
        if ( $line =~ /^\>([\w:-])+/ ) {
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

1;

__END__

=head1 SYNOPSIS

    use AlignDB::Util qw(:all);

