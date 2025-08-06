# SLURM/SBATCH Error Handling & Debugging 

## Dashing through the -e
Using `#!/bin/bash -e` (or `set -e` within the script) tells Bash to exit immediately if any command returns a non-zero (i.e., failure) status - **and can also It also serve as a mechanism for ensuring your failed jobs show a status of `FAILED` in `sacct` output. Here’s a breakdown of the pros and cons of enabling that behavior in your BASH/SLURM submission scripts.

### Benefits of -e
- **Fast failure detection**  
  Ensures your script halts as soon as something goes wrong, preventing it from blindly continuing and possibly producing unreliable output.  
- **Easier debugging**  
  By stopping at the first failing command, your error logs point straight to the root problem rather than burying it amidst downstream noise.  
- **Safer automation**  
  In automated pipelines, it's often better to fail fast than to carry on with partial or corrupted results—especially in large-scale data processing.  
- **Reduces masking of errors**  
  Comparatively, scripts without -e might silently skip failures and mislead users by generating misleading success indications.  
- **Clear Separation of stdout and stderr**  
  Keeping error logs separate from standard output allows you to easily spot issues (e.g., exceptions or warnings) without sifting through normal messages.  
  ```bash
  #SBATCH -o job_%j.out
  #SBATCH -e job_%j.err
  ```
- **Cleaner and More Focused Debugging**  
  When a job fails, you can go directly to the .err file, which contains only error messages and stack traces—ideal for quickly diagnosing problems.  
- **Organizational Clarity**  
  Using distinct filenames like align_12345.err and align_12345.out helps you associate each file with the corresponding job and purpose, making logs easy to manage and store.  
- **Support for Automation**  
  Scripts can programmatically parse .err files to send alerts or trigger follow-up tasks only when errors are detected, rather than parsing combined logs.  

### Drawbacks of -e
- **Too aggressive for some commands** 
  Commands expected to fail (like grep searching non-existent matches) or transient network operations may prematurely stop entire pipelines.  
- **Lack of granular control**  
  Unless you explicitly guard individual commands (e.g., using || true or checking if ...; then), the script won't continue past any error—even if it's non-critical.
- **Silent termination**  
  The script terminates immediately, often without obvious messaging unless you've instrumented it with echo statements or logging steps before each critical command.  
- **Increased File Clutter**  
  You'll generate at least two files per job—stdout and stderr—which can clutter directories, especially with many small jobs.  
- **Mixed Message Ordering**  
  Errors and log messages might occur out of real-time order when split between two files. This could confuse debugging if a failure depends on interleaving logic.  

*In situations where you’re not sure . . .*  
- **Use in combination with control structures** - Rather than a blanket `-e`, wrap risky commands in conditional logic ( `||` = “or” ):
  ```bash
  #!/bin/bash -e
  
  echo "Step 1: critical task"
  /path/to/critical_command
  
  echo "Step 2: optional task"
  /path/to/optional_command || echo "Warning: optional task failed, continuing . . ."
   ```
  
- **Alternatively**: Set `-e` later:
  ```bash
  #!/bin/bash
  # Preliminary setup steps not strict
  module load module1
  module load module2
  
  # From here on, fail fast
  set -e
  
  critical_command1
  critical_command2
  ```
  
- **Use `trap`** to improve error reporting:
  ```bash
  #!/bin/bash -e
  trap 'echo "Error on line $LINENO"; exit 1' ERR
  
  echo "Running step 1"
  do_something
  
  echo "Running step 2"
  do_something_else
  ```
  
### When to use for your HPC workflows:
- **Use it** in production or data-processing jobs where any mistake should abort execution immediately.  
- **Avoid it** (or guard commands) when scripts contain:  
  - Inner conditionals  
  - Commands that may fail and return non-zero on non-critical issues  
  - Loops or `grep | wc`  pipelines expected to produce empty results  
               
### Alternatives:
- **Naming Patterns**: Use `%x` (job name) and `%j` (job ID) in filenames:
  ```bash
  #SBATCH -o %x_%j.out
  #SBATCH -e %x_%j.err
  ```
- **Conditional Error Handling**: In your script, check if the stderr file has content and, if so, trigger notifications:
  ```bash
  if [[ -s "${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err" ]]; then
    echo "Job failed—see STDERR log" | mail -s "Job error" first.last@domain.com
  fi
  ```
- **Cleanup Strategy**: Automatically compress or delete old log files to manage filesystem usage.  
- **Merged Output Only**: Skip `-e`, so both `stdout` and `stderr` go to the same file via `--output=job_%j.log`.  

---

## Requeue’in the Night Away
When you submit a job with the `--requeue` flag (e.g., `#SBATCH --requeue`), you're telling SLURM: if this job is terminated prematurely—such as due to node failure, preemption, or reaching its time limit—it should be automatically requeued rather than being marked as failed.

### Benefits of --requeue
- **Resilience to interruptions**: If a node crashes or your job is preempted (especially on preemptible partitions), SLURM will put it back in the queue with the same job ID and retained priority, allowing it to run again cleanly. No need to re-submit manually.  
- **Preserves job accounting & priority**: The requeued job uses the same job ID, making resource tracking consistent. It also remains in your place in the scheduling queue, avoiding delays caused by losing your priority position.  
- **Supports checkpoint-restart workflows**: Combined with robust checkpointing (via software or custom scripts triggered by `SIGTERM`, usually around 10 seconds before timeout), jobs can resume from intermediate state rather than starting over.  

### Drawbacks of --requeue
- **Requires proper checkpointing**: Without saving progress before termination, jobs will restart from scratch—wasting time and resources. You must integrate checkpoint logic into your application or scripts.  
- **May lead to endless loops**: If the requeued job hits the same error repeatedly (e.g., due to bad input or configuration), it could cycle forever—unless you include safeguards or limits on retries.  
- **Resource availability and queue delays**: A requeued job re-enters the queue, so it might wait behind other jobs. For time-sensitive tasks, this can hurt overall throughput.  
- **Complex error diagnosing**: Since job restarts can trigger without your intervention, logs may span multiple runs. You’ll need to check logs from each attempt to debug issues, which can be more complex.  

### Best Practices with --requeue
- **Checkpoint smart**: Use `--signal=B:SIGTERM@60` or similar to get a `SIGTERM` some time before SLURM kills your job. In a trap handler, save intermediate state so you can restart from the last checkpoint.  
- **Set requeue limits** (`MAX=5` on our clusters): Either via in-script counters or via SLURM configuration, prevent infinite restarts by limiting retries.  
- **Combine with arrays/dependencies**: If you’re running arrays, you can requeue only failed elements. Or structure workflows so each step requeues on failure, separate from others.  
- **Log clearly**: Include timestamps and a “run attempt number” in your logs. For example: `echo "RUN #$SLURM_RESTART_COUNT at $(date)"`  
- **Test interactively**: Before enabling auto-requeueing, test checkpoint-restart logic in an interactive session (e.g., `srun --requeue`), so you validate correct behavior.  

#### Use --requeue when:
- You expect possible preemption or node/network failure.
- You have checkpointing built into your workflow or application.
- You want cleaner accounting and job management.

#### Avoid or use with caution when:
- Jobs are short or not checkpointable.
- You risk creating infinite restart loops.
- Scheduling delays would harm your workflow.

---

## Snappin' & Trappin'
Here are a few polished examples showing how to use `trap` statements in `SBATCH` scripts to improve error reporting and handling — especially helpful for bioinformatics pipelines running on SLURM. Using well-crafted trap statements in SBATCH scripts boosts robustness by:  
- Providing useful error context automatically.  
- Preventing data loss by checkpointing when jobs are interrupted.  
- Enabling automated and successful job recovery with SLURM’s requeue support.  

Included are two real-world examples that combine **trap statements**, **checkpointing**, and SLURM-**requeue** support in *Snakemake* and *Nextflow* workflows. These scripts are also a part of our “user utilities” on both clusters: `/data/user-utils/user-scripts/`

### Trap ERR for Immediate Error Identification
To catch any unexpected error and report exactly where it occurred, combine `set -eEuo pipefail` with an `ERR` trap:

  `errtrap.sbatch`
  ```bash
  #!/bin/bash
  #SBATCH --job-name=issa-trap
  #SBATCH --output=issa-trap_%j.out
  #SBATCH --error=issa-trap_%j.err
  
  set -Eeuo pipefail
  trap 'echo "[$(date)] ERROR in ${BASH_SOURCE[0]} at line ${BASH_LINENO[0]}: command \"${BASH_COMMAND}\" failed."; exit 1' ERR
  
  # Load necessary modules
  ml bwa
  ml samtools
  
  # Create dummy reference genome
  cat <<EOF > ref.fa
  >chr1
  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  EOF
  
  # Create dummy FASTQ reads
  cat <<EOF > reads.fq
  @read1
  ACGTACGTACGT
  +
  FFFFFFFFFFFF
  @read2
  ACGTACGTACGT
  +
  FFFFFFFFFFFF
  EOF
  
  # Build BWA index (required for bwa mem)
  bwa index ref.fa
  
  # Run alignment and sort the result
  bwa mem ref.fa reads.fq > aln.sam
  samtools sort -O BAM aln.sam -o aln.sorted.bam
  ```
  ```bash
  $ sbatch issa-trap_25218
  Submitted batch job 25218
  ```
  ```bash
  $ cat issa-trap_252182.err
  [bwa_index] Pack FASTA... 0.00 sec
  [bwa_index] Construct BWT for the packed sequence...
  [bwa_index] 0.00 seconds elapse.
  [bwa_index] Update BWT... 0.00 sec
  [bwa_index] Pack forward-only FASTA... 0.00 sec
  [bwa_index] Construct SA from BWT and Occ... 0.00 sec
  [main] Version: 0.7.17-r1188
  [main] CMD: bwa index ref.fa
  [main] Real time: 0.008 sec; CPU: 0.002 sec
  [M::bwa_idx_load_from_disk] read 0 ALT contigs
  [M::process] read 2 sequences (24 bp)...
  [M::mem_process_seqs] Processed 2 reads in 0.000 CPU sec, 0.000 real sec
  [main] Version: 0.7.17-r1188
  [main] CMD: bwa mem ref.fa reads.fq
  [main] Real time: 0.002 sec; CPU: 0.002 sec
  ```

**How it works**:
- `set -Eeuo pipefail`:  
  - `-e`: exit on error  
  - `-E`: propagate ERR traps into subshells/functions  
  - `-u`: error on unset variables  
  - `pipefail`: fail on pipeline error  
- The `trap` prints the script name, line number (top of the faulty frame), and the failed command—then exits.  
- Clear, automatic failure reporting with context—no guessing which command failed.  

### Trap SIGTERM for Graceful Checkpointing
If your job uses `--requeue`, SLURM may send `SIGTERM` before requeuing. Trap it to save checkpoints before exit:

  `graceful-checkpoint.sbatch`
  ```bash
  #!/bin/bash
  #SBATCH --job-name=graceful-checkpoint
  #SBATCH --output=graceful-ckpt_%j.out
  #SBATCH --error=graceful-ckpt_%j.err
  #SBATCH --time=2:00:00
  #SBATCH --signal=B:USR1@60        # send SIGUSR1 60s before timeout
  #SBATCH --requeue
  
  set -eEuo pipefail
  trap 'echo "CATCHED SIGUSR1; saving checkpoint..."; save_checkpoint; exit 1' USR1
  
  ml snakemake
  
  save_checkpoint() {
    echo "Snakemake checkpoint at $(date)" >> checkpoint_${SLURM_JOB_ID}.log
    cp -r .snakemake .snakemake_checkpoint_${SLURM_JOB_ID}
  }
  
  # Resume if checkpoint exists
  if [[ -f checkpoint_${SLURM_JOB_ID}.dat ]]; then
    echo "Resuming from checkpoint"
    # resume logic...
  else
    echo "Starting fresh"
  fi
  
  # Example workload
  for i in $(seq 1 100); do
    echo "Step $i"
    sleep 30
    echo $i > checkpoint_${SLURM_JOB_ID}.dat
  done
  
  echo "Completed successfully"
  rm -f checkpoint_${SLURM_JOB_ID}.dat
  ```
  ```bash               
  $ sbatch graceful-checkpoint.sbatch
  Submitted batch job 252184
  ```
  ```bash
  $ cat graceful-ckpt_252184.out
  Starting fresh
  Step 1
  Step 2
  Step 3
  ```

**How it works**:
- `--signal=B:USR1@60` sends `SIGUSR1` 60 seconds before reaching the time limit.
- The `trap` handles the signal and runs save_checkpoint so work isn't lost.
- `--requeue` ensures SLURM re-submits the job.
- On restart, the script checks for the checkpoint and resumes rather than starting over.

#### Best Practices:
- Use `-Ee` to **ensure `ERR` traps fire** even in functions/subshells: `set -Eeuo pipefail`
- **Print informative context**: Use `${BASH_SOURCE[0]}`, `${BASH_LINENO[0]}`, and `${BASH_COMMAND} `for precise debugging output.
- **Checkpoint on signals**: Common choices include USR1 (before time limit) and TERM (if node failure). Use:
trap 'save_checkpoint; exit 1' USR1 TERM
- **Combine with --requeue**: Ensure SLURM can requeue your job after trapping the signal and checkpointing.

### Snakemake + SBATCH Script with Trap & Requeue
Use a wrapper SBATCH script to invoke Snakemake and capture signals (you will need a valid `Snakefile` for this script to run):

  `snakemake-on-a-plane.sbatch`
  ```bash
  #!/bin/bash
  #SBATCH --job-name=snakemake-on-a-plane
  #SBATCH --output=snake_wrap_%j.out
  #SBATCH --error=snake_wrap_%j.err
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=4
  #SBATCH --mem=16G
  #SBATCH --time=24:00:00
  #SBATCH --requeue
  #SBATCH --signal=B:USR1@300
  
  set -Eeuo pipefail
  trap 'echo "[$(date)] Signal USR1 received—checkpointing and requeueing..."; snakemake --rerun-incomplete; exit 1' USR1 ERR
  
  ml snakemake
  
  # Run Snakemake with SLURM integration
  snakemake \
    --snakefile Snakefile \
    --cores 4 \
    --jobs 100 \
  #  --cluster "sbatch --parsable --partition=all --cpus-per-task={threads} --mem={resources.mem_mb//1024}G --time={resources.runtime_min//60}:00" \
  #  --latency-wait 60 \
  #  --rerun-incomplete
  ```
  ```
  Snakefile
  # Snakefile
  
  rule all:
      input:
          "results/hello.txt"
  
  rule say_hello:
      output:
          "results/hello.txt"
      shell:
          "echo 'GET THESE SNAKES OFF THIS PLANE!' > {output}"
  ```
  ```bash               
  $ sbatch snakemake-on-a-plane.sbatch
  Submitted batch job 252191
  ```
  ```bash
  $ cat results/hello.txt
  GET THESE SNAKES OFF THIS PLANE!
  ```
  ```bash
  $ cat snake_wrap_252191.err
  Assuming unrestricted shared filesystem usage.
  Building DAG of jobs...
  Using shell: /usr/bin/bash
  Provided cores: 4
  Rules claiming more threads will be scaled down.
  Job stats:
  job          count
  --------- ** -------
  all              1
  say_hello        1
  total            2
  
  Select jobs to execute...
  Execute 1 jobs...
  
  [Fri Jul 18 12:37:38 2025]
  localrule say_hello:
      output: results/hello.txt
      jobid: 1
      reason: Missing output files: results/hello.txt
      resources: tmpdir=/path/to/host
  
  [Fri Jul 18 12:37:38 2025]
  Finished job 1.
  1 of 2 steps (50%) done
  Select jobs to execute...
  Execute 1 jobs...
  
  [Fri Jul 18 12:37:38 2025]
  localrule all:
      input: results/hello.txt
      jobid: 0
      reason: Input files updated by another job: results/hello.txt
      resources: tmpdir=/path/to/host
  
  [Fri Jul 18 12:37:38 2025]
  Finished job 0.
  2 of 2 steps (100%) done
  Complete log: .snakemake/log/2025-07-18T123738.691055.snakemake.log
  ```

**How it works**:
- Captures `SIGUSR1` (sent 5 min before timeout) and re-runs incomplete tasks via Snakemake’s `--rerun-incomplete` mechanism.
- Uses `--requeue`, retaining the same job ID and priority if interrupted.
- Snakemake's `--cluster` mode automatically submits individual rules as SLURM jobs, mapping `threads`, `mem_mb`, and `runtime_min` to `SBATCH` directives

### Nextflow + Checkpoint Trap + srun
Wrap a Nextflow pipeline in `SBATCH`, enabling graceful checkpoint handling when the job is signaled (you’ll need a valid `main.nf` to run this script):

  `nextflowing.sbatch`
   ```bash
  #!/bin/bash
  #SBATCH --job-name=nextflowing
  #SBATCH --output=nextflow_wrap_%j.out
  #SBATCH --error=nextflow_wrap_%j.err
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=4
  #SBATCH --mem=16G
  #SBATCH --time=48:00:00
  #SBATCH --requeue
  #SBATCH --signal=B:USR1@300
  
  set -Eeuo pipefail
  
  cancel_nextflow() {
    echo "[$(date)] Caught SIGUSR1—terminating Nextflow gracefully"
    srun --jobid="${SLURM_JOB_ID}.0" kill -TERM "$NEXTFLOW_PID"
  }
  trap cancel_nextflow USR1 ERR
  
  ml nextflow
  
  # Launch Nextflow as a background job using srun
  srun --exclusive nextflow run main.nf -profile slurm &
  NEXTFLOW_PID=$!
  
  wait "$NEXTFLOW_PID"
  exit_status=$?
  
  if [[ $exit_status -ne 0 ]]; then
    echo "[$(date)] Nextflow exited with $exit_status—requeueing..."
    exit 1
  fi
  
  echo "[$(date)] Pipeline completed successfully. If you are quiet enough, you will hear the flow of the universe."
  exit 0
  ```
  ```
  main.nf
  nextflow.enable.dsl=2
  
  process sayHi {
      output:
      file 'hello.txt'
  
      script:
      '''
      echo "Hi from Nextflow on SLURM!" > hello.txt
      '''
  }
  
  workflow {
      sayHi()
  }     
  ```
  ```bash
  $ sbatch nextflowing.sbatch
  Submitted batch job 252200
  ```
  ```bash
  $ cat nextflow_wrap_252200.err
  Nextflow 25.04.6 is available - **Please consider updating your version to it
  ```
  ```bash
  $ cat nextflow_wrap_252200.out
  
  N E X T F L O W   ~  version 24.10.0
  
  Launching `main.nf` [berserk_brenner] DSL2 - **revision: d21dfb9617
  
  [- **       ] sayHi -
  
  executor >  local (1)
  [0c/03b0a2] sayHi | 0 of 1
  
  executor >  local (1)
  [0c/03b0a2] sayHi | 1 of 1 ✔
  
  [Fri Jul 18 03:31:22 PM MDT 2025] Pipeline completed successfully. If you are quiet enough, you will hear the flow of the universe.
  ```

Key features:
- Runs Nextflow under `srun` so it resides in a SLURM job step and can receive signals correctly.
- Captures `SIGUSR1`, sends `TERM` to the Nextflow process to allow graceful shutdown and checkpoint cleanup.
- With `--requeue`, the script will be resubmitted with the same job ID if terminated—enabling Nextflow’s `-resume` capability from its last checkpoint.

Snakemake handles clusters natively by creating per-rule `SBATCH` jobs. Wrapping it allows graceful response to resource signals and safe requeueing. Nextflow doesn’t auto-handle SLURM signals, but using an SBATCH wrapper with srun lets you catch signals and terminate Nextflow for checkpointing and restart. Both approaches maintain reproducibility while leveraging SLURM's `--signal` + `--requeue` mechanism for robust, long-running workloads.
