source /rds/user/fd385/hpc-work/software/mace-venv/bin/activate
path=/rds/user/fd385/hpc-work/software/mace/mace/cli

python  $path/run_train.py \
    --name="MACE_model" \
    --foundation_model="large" \
    --default_dtype='float32'\
    --r_max=6.0 \
    --valid_fraction=0.1 \
    --E0s="{1:-0.14003657,8:-0.01132296}" \
    --loss="stress" \
    --train_file="training_dataset.xyz" \
    --scaling="rms_forces_scaling"\
    --error_table="PerAtomMAE"\
    --batch_size=2 \
    --max_num_epochs=2000 \
    --forces_weight=550 \
    --energy_weight=10000 \
    --swa \
    --start_swa=1000 \
    --swa_forces_weight=550 \
    --swa_energy_weight=10000 \
    --swa_stress_weight=5 \
    --swa_lr=0.00025 \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --restart_latest \
    --device=cuda \
    --seed=1 \
