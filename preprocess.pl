#!/usr/bin/env perl

############################################################################
## Preprocesses an input graph of undirected links to a format required by
## graph statistics programs; namely,
##
##  1) each undirected link is converted to a pair of directed links, and
##  2) links are grouped by the source node.
##
## Property 2 ensures that all outgoing links for a given node appear
## in consecutive lines in the output graph file.  We don't actually
## have to sort the output lines to achieve property 2, but we do so
## anyway in case it proves useful.
##
## ---------------------------------------------------------------------
## Copyright (C) 2010 The Regents of the University of California.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

open SORT, "| sort -k 1,2n -k 2,3n"
  or die "ERROR: couldn't open pipe to 'sort' command";

while (<>) {
  next if /^\s*$/ || /^\#/;

  if (/^(\d+)\s+(\d+)\s*$/) {
    my $s = $1;
    my $d = $2;
    if ($s == $d) {
      warn "discarding self loop on node $s\n";
      next;
    }

    print SORT "$s $d\n$d $s\n";
  }
  else {
    die "ERROR: malformed line: $_";
  }
}

