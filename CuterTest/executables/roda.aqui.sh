#!/bin/bash

# EDIT dcicpp.spc manually

MAXIT='200000'
MAXTIME='600'

rm -f latex*
./runlist.sh TestLists/testes.equ > testes.equ.out.$MAXIT.$MAXTIME
wc -l latex* > testes.equ.wc.$MAXIT.$MAXTIME
rm -f latex*
./runlist.sh TestLists/testes.ineq > testes.ineq.out.$MAXIT.$MAXTIME
wc -l latex* > testes.ineq.wc.$MAXIT.$MAXTIME
rm -f latex*
./runlist.sh TestLists/testes.unc > testes.unc.out.$MAXIT.$MAXTIME
wc -l latex* > testes.unc.wc.$MAXIT.$MAXTIME
rm -f latex*
./runlist.sh TestLists/testes.gencon > testes.gencon.out.$MAXIT.$MAXTIME
wc -l latex* > testes.gencon.wc.$MAXIT.$MAXTIME

