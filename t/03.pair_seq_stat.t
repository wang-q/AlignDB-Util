#!/usr/bin/perl
use strict;
use warnings;

use Test::More qw(no_plan);
use Test::Number::Delta within => 1e-2;

BEGIN {
    use_ok( 'AlignDB::Util', qw(pair_seq_stat multi_seq_stat) );
}

{
    print "#pair_seq_stat\n";

    #$seq_legnth,            $number_of_comparable_bases,
    #$number_of_identities,  $number_of_differences,
    #$number_of_gaps,        $number_of_n,
    #$number_of_align_error, $pi,
    #$first_seq_gc,          $average_gc,
    my @data = (

        #AAAATTTTGG
        #AAAATTTTTG
        [   [qw{ AAAATTTTGG AAAATTTTTG }],
            [ 10, 10, 9, 1, 0, 0, 0, 0.1, 0.2, 0.15 ],
        ],

        #TTAGCCGCTGAGAAGC
        #GTAGCCGCTGA-AGGC
        [   [qw{ TTAGCCGCTGAGAAGC GTAGCCGCTGA-AGGC }],
            [ 16, 15, 13, 2, 1, 0, 0, 0.1333, 0.5625, 0.6146 ],
        ],

        #GATTATCATCACCCCAGCCACATA
        #GATTTT--TCACTCCATTCGCATA
        [   [qw{ GATTATCATCACCCCAGCCACATW GATTTT--TCACTCCATTCGCATA }],
            [ 24, 21, 16, 5, 2, 1, 0, 0.2381, 0.4783, 0.4209 ],
        ],

    );

    for my $i ( 0 .. @data - 1 ) {
        my ( $seq_pair_ref, $except_ref ) = @{ $data[$i] };
        my $result_ref = pair_seq_stat(@$seq_pair_ref);
        print "\n";
        delta_ok( $result_ref, $except_ref, "stat $i" );

        my $result_multi_ref = multi_seq_stat(@$seq_pair_ref);
        print "\n";
        delta_ok( $result_multi_ref, $except_ref, "stat $i" );
    }
}

