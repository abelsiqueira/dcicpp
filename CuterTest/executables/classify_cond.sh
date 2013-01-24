#!/bin/bash
for type in assert conv maxrest maxiter rhomax unlimited restfail short time infeasible
do
  grep $type cond.inf > cond_$type.inf
done
wc -l cond_*
