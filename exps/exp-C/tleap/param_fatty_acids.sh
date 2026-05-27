#!/bin/bash
# Parameterize fatty acid chains with GAFF2 + AM1-BCC
set -e
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap

echo "=== Parameterizing C16 monoacid (palmitate) ==="
echo "CCCCCCCCCCCCCCCC(=O)[O-]" > /tmp/fa.smi
antechamber -i /tmp/fa.smi -fi smi -o c16_monoacid.mol2 -fo mol2 -c bcc -at gaff2 -nc -1 -pf yes 2>&1 | grep -E "(Info|Error|Fatal|charge)" | head -10
echo "C16 monoacid done."

echo "=== Parameterizing C18 monoacid (stearate) ==="
echo "CCCCCCCCCCCCCCCCCC(=O)[O-]" > /tmp/fa.smi
antechamber -i /tmp/fa.smi -fi smi -o c18_monoacid.mol2 -fo mol2 -c bcc -at gaff2 -nc -1 -pf yes 2>&1 | grep -E "(Info|Error|Fatal|charge)" | head -10
echo "C18 monoacid done."

echo "=== Parameterizing C18 diacid (octadecanedioate, -2) ==="
echo "O=C([O-])CCCCCCCCCCCCCCCCC(=O)[O-]" > /tmp/fa.smi
antechamber -i /tmp/fa.smi -fi smi -o c18_diacid.mol2 -fo mol2 -c bcc -at gaff2 -nc -2 -pf yes 2>&1 | grep -E "(Info|Error|Fatal|charge)" | head -10
echo "C18 diacid done."

echo "=== Parameterizing C16 diacid (hexadecanedioate, -2) ==="
echo "O=C([O-])CCCCCCCCCCCCCCC(=O)[O-]" > /tmp/fa.smi
antechamber -i /tmp/fa.smi -fi smi -o c16_diacid.mol2 -fo mol2 -c bcc -at gaff2 -nc -2 -pf yes 2>&1 | grep -E "(Info|Error|Fatal|charge)" | head -10
echo "C16 diacid done."

ls -lh c16_monoacid.mol2 c18_monoacid.mol2 c18_diacid.mol2 c16_diacid.mol2 2>/dev/null
echo "Parameterization complete."
