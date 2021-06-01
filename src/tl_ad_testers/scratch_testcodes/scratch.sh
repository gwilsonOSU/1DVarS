#!/bin/zsh

fin=tl_qtrans_vanderA.m
grep '^tl_' $fin | \
    sed 's/^\(.*\)=.*/\1/g'
