#!/usr/bin/perl
use strict;
use warnings;

use Test::More qw(no_plan);

BEGIN {
    use_ok( 'AlignDB::Util', qw(calc_gc_ratio) );
}

{
    print "#seq_length\n";

    my @seqs = (
        [ "ATAA",            0 ],
        [ "AtaA",            0 ],
        [ "CCGC",            1 ],
        [ "CcGc",            1 ],
        [ "TAGggATaaC",      0.4 ],
        [ "GCaN--NN--NNNaC", 0.6 ],
    );

    for my $i ( 0 .. @seqs - 1 ) {
        my ( $ori, $expected ) = @{ $seqs[$i] };
        my $result = calc_gc_ratio($ori);
        print "\n";
        print "original: $ori\n";
        is( $result, $expected, "calc_gc_ratio_$i" );
    }
}

