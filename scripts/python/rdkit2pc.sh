# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/research/users/amirhs/miniconda3_2/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
eval "$__conda_setup"
else
  if [ -f "/research/users/amirhs/miniconda3_2/etc/profile.d/conda.sh" ]; then
. "/research/users/amirhs/miniconda3_2/etc/profile.d/conda.sh"
else
  export PATH="/research/users/amirhs/miniconda3_2/bin:$PATH"
fi
fi
unset __conda_setup
conda activate my-rdkit-env
python names2smiles.py cas.txt smiles.txt
python pc_rdkit.py
