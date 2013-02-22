#!/bin/bash

conv=$(wc -l latex_* ../latex_* | grep conv | awk '{s=s+$1}END{print s}')
total=$(wc -l latex_* ../latex_* | grep total | awk '{s=s+$1}END{print s}')

echo -n "$conv / $total = "
echo "scale=5; $conv/$total" | bc
