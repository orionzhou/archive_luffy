#!/usr/bin/perl -w
use strict;
use DBI;
our $dbh = DBI->connect("DBI:mysql:database=mt;host=202.117.161.53", "genius", "prodigy", {"RaiseError"=>1});