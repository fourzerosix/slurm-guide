## Introduction to SBATCH
As a workload manager, [SLURM](https://slurm.schedmd.com/documentation.html) allocates compute resources (nodes, CPUs, GPUs, memory) and provides a framework for starting and monitoring jobs on those resources. Users submit batch jobs via [SBATCH](slurm.schedmd.com/sbatch.html) scripts, which are ordinary BASH scripts that include special #SBATCH directives (comments) telling SLURM what resources and options to use. A basic SBATCH script looks like:
  ```bash
  #!/bin/bash
  #SBATCH --job-name=my_job       # Job name
  #SBATCH --output=out_%j.txt     # Standard output (%j expands to job ID)
  #SBATCH --error=err_%j.txt      # Standard error
  #SBATCH --partition=all         # Partition (queue) name
  #SBATCH --ntasks=1              # Number of tasks (MPI ranks) 
  #SBATCH --cpus-per-task=1       # Number of CPU cores per task
  #SBATCH --mem=1G                # Total memory
  #SBATCH --time=00:10:00         # Walltime (D-HH:MM:SS)
  
  # (Optional) load modules or activate envs here
  # module load mytool
  
  # Commands for your job:
  echo "Running on $(hostname)"
  ```

In this template, each `#SBATCH` line specifies an option: e.g. `--job-name` names the job, `--output`/`--error` set log files, `--partition` selects a queue, and `--time` sets the max run time. These directives are simply comments in the script that SLURM reads to configure the job. After the directives, the script runs normal shell commands (e.g. loading software modules, activating Conda environments, then executing code). To submit such a script, run `sbatch script_name.[sh:sbatch]`; SLURM will queue and then execute it when resources are available, returning a job ID.

SLURM provides various commands to monitor and control jobs. Commonly used tools include `sinfo` (to view node/partition status), `squeue` (to see queued/running jobs), `scancel` (to kill jobs), and `sacct` (to report on completed jobs). For example, `squeue -u $USER` shows your jobs, and `sacct -j JOBID` shows history and resource usage for a finished job. Together, SBATCH scripts and SLURM commands form the basis of batch processing on the cluster.

> Interactive Session
> Sometimes you need an interactive shell on a compute node (for testing or GUI apps). You can allocate resources with srun or salloc. For example:
  ```bash
  srun --partition=all --ntasks=1 --cpus-per-task=2 --mem=4G --time=01:00:00 --pty bash
  ```

This command requests an interactive shell ( `--pty bash` ) on a compute node for 1 hour with 2 CPUs and 4 GB memory. Once granted, you get a shell prompt on a node and can run commands as if logged in. (On some systems you may use `salloc` similarly.) In many clusters, there may also be a helper script (like interact) that wraps srun. Interactive sessions are useful for debugging jobs or running GUI programs on the node.

---

### Example SBATCH Job Scripts
Below are annotated SBATCH script examples for common use cases. Each script begins with `#!/bin/bash` and `#SBATCH` directives; adjust partition names or resource amounts for your cluster (“`all`”, “`gpu`”, “`himem`” partitions; *max walltime is **28** days*).

#### Standard Single-Core CPU Job
This simplest job uses one CPU on the default “ `all`” partition. It runs a single command ( `echo` ) and writes logs to files.
  ```bash
  #!/bin/bash
  #SBATCH --job-name=!NT3R3$T!NG_T1M35        # Job name
  #SBATCH --output=interesting-times_%j.out   # STDOUT (%j = JobID)
  #SBATCH --error=interesting-times_%j.err    # STDERR
  #SBATCH --partition=all                     # Partition (all = default)
  #SBATCH --ntasks=1                          # Total tasks (MPI ranks)
  #SBATCH --cpus-per-task=1                   # CPU cores per task
  #SBATCH --mem=1G                            # Total memory (1 GB)
  #SBATCH --time=00:01:00                     # Time limit (1 minute)

  echo "May you live in interesting times. Regards, $(hostname)"
  ```
  ```bash
  $ sbatch single-core-cpu.sbatch
  Submitted batch job 237762
  ```
  ```bash
  $ cat interesting-times_237762.out
  May you live in interesting times. Regards, host
  ```

**Explanation**: *This script requests 1 node, 1 task, and 1 CPU core ( --ntasks=1, --cpus-per-task=1 ) with 1 GB RAM and 1 min walltime. The echo command will run on a single core. SLURM will create output/error files ( here named interesting-times _<jobid>.[out:err] ) where all output is logged. The directives used ( --job-name, --output, --partition, --ntasks, --cpus-per-task, --time ) are standard SBATCH options. When submitted, SLURM prints a job ID; use squeue or scontrol show job <ID> to track it.*

---

#### Multi-Node MPI Job
For distributed-memory [MPI](https://www.open-mpi.org/) jobs, request multiple nodes and tasks. For example, an MPI job using 4 processes across 2 nodes:
  ```bash
  #!/bin/bash
  #SBATCH --job-name=P@$$!NG_M3SS@G3S
  #SBATCH --output=passing-messages_%j.out
  #SBATCH --partition=all
  #SBATCH --nodes=2
  #SBATCH --ntasks=2
  #SBATCH --ntasks-per-node=1
  #SBATCH --cpus-per-task=4      # 4 CPU cores for threading
  #SBATCH --mem=4G               # (Optional - have encountered issues without a --mem= specification when using MPI
  #SBATCH --time=00:01:00
  
  ml openmpi/5.0.7-rwefebl
  ml osu-micro-benchmarks/7.5-23m3nxv
  
  srun --mpi=pmix osu_bw
  ```
  ```bash
  $ sbatch mmmpi.sbatch
  Submitted batch job 238040
  ```
  ```bash
  $ squeue --me --format="%8i %14j %22S %10M %6D %6C %12m %10T %R"
  JOBID    NAME           START_TIME             TIME       NODES  CPUS   MIN_MEMORY   STATE      NODELIST(REASON)
  238040   P@$$!NG_M3SS@G 2025-06-27T12:52:28    0:02       2      8      4G           RUNNING    ai-rmlc[08-09]
  ```
  ```bash
  $ cat passing-messages_238040.out
  # OSU MPI Bandwidth Test v7.5
  # Datatype: MPI_CHAR.
  # Size      Bandwidth (MB/s)
  1                       1.97
  2                       3.95
  4                       7.93
  8                      15.89
  16                     31.66
  32                     63.09
  64                    121.81
  128                   245.92
  256                   467.41
  512                   916.54
  1024                 1785.72
  2048                 3597.19
  4096                 5928.28
  8192                 8834.14
  16384                9627.58
  32768                9931.16
  65536               13374.24
  131072              15678.86
  262144              16688.54
  524288              16395.74
  1048576             17055.63
  2097152             16900.65
  4194304             17873.02
  ```

**Explanation**: *This requests 2 nodes ( `--nodes=2` ) and a total of 4 MPI tasks  ( `--ntasks=2` ); by default, SLURM will distribute the 2 tasks across the 2 nodes. We set `--cpus-per-task=4` for four CPU-cores/threading. After loading MPI modules, we run the program with `srun` (or `mpirun`); SLURM then starts 2 processes across the allocated nodes. Here we ask for 4 CPUs on each node ( `--cpus-per-task=4`, `--ntasks=2` means two tasks that can each use 4 threads ). We load any needed modules (e.g., `openmpi`/`mpi` application). The example program is launched with `srun`, which distributes it onto the allocated cores.*

---

#### GPU Job
For GPU-accelerated jobs, request GPU resources on the GPU partition. For example, to run on 1 GPU:
  ```bash
  #!/bin/bash
  #SBATCH --job-name=i_like_pictures
  #SBATCH --output=pretty-pictures_%j.out
  #SBATCH --partition=gpu      # GPU partition (if “gpu” is the correct name)
  #SBATCH --nodes=1
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=2
  #SBATCH --gres=gpu:1         # Request 1 GPU
  #SBATCH --mem=16G
  #SBATCH --time=01:00:00

  ml cuda       # load CUDA or other GPU driver modules as needed

  srun nvidia-smi
  ```
  ```bash
  $ sbatch gpu-job.sbatch
  Submitted batch job 238098
  ```
  ```bash
  $ cat pretty-pictures_238098.out
  Fri Jun 27 13:13:24 2025
  +---------------------------------------------------------------------------------------+
  | NVIDIA-SMI 545.23.08              Driver Version: 545.23.08    CUDA Version: 12.3     |
  |-----------------------------------------+----------------------+----------------------+
  | GPU  Name                 Persistence-M | Bus-Id        Disp.A | Volatile Uncorr. ECC |
  | Fan  Temp   Perf          Pwr:Usage/Cap |         Memory-Usage | GPU-Util  Compute M. |
  |                                         |                      |               MIG M. |
  |=========================================+======================+======================|
  |   0  NVIDIA A100-PCIE-40GB          On  | 00000000:E2:00.0 Off |                    0 |
  | N/A   29C    P0              33W / 250W |      4MiB / 40960MiB |      0%      Default |
  |                                         |                      |             Disabled |
  +-----------------------------------------+----------------------+----------------------+
  ```

**Explanation**: *We specify `--partition=gpu` to use GPU nodes. `--gres=gpu:1` asks SLURM to allocate one GPU on the node. We also allocate a couple of CPU cores and some RAM. We then load any GPU drivers/modules (for example, CUDA) and run the program that uses the GPU.*

#### High-Memory Job (with Miniconda3)
For tasks needing a large amount of RAM, use the “`himem`” partition and request more memory. This job loads miniconda3, creates a [Matplotlib](https://matplotlib.org/) conda environment from our in-house conda repo ( `/path/to/conda-envs/` ) - and runs a numpy large matrix operation (this will produce a new conda environment for you, along with a histogram and log file):
  ```bash
  #!/bin/bash
  #SBATCH --job-name=matplotlib_himem
  #SBATCH --output=matplotlib_%j.out
  #SBATCH --error=matplotlib_%j.err
  #SBATCH --partition=himem
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=4
  #SBATCH --mem=256G
  #SBATCH --time=24:00:00
  
  # robust bash behavior
  set -Eeuo pipefail
  trap 'echo "[$(date)] Error in ${BASH_SOURCE[0]} at line ${BASH_LINENO[0]}: \"$BASH_COMMAND\" failed." >&2; exit 1' ERR
  
  echo "[$(date)] Starting job on $(hostname). SLURM_JOB_ID=$SLURM_JOB_ID"
  
  # load Miniconda3
  module load miniconda3        # adjust to exact module name/version if needed
  
  # create env only if not already existing and sentient
  ENV_NAME=matplotlib_3.10.0
  YAML_FILE=/path/to/conda-envs/matplotlib_3.10.0.yaml
  
  # set up env path in project/scratch directory
  ENV_DIR=$HOME/conda_envs/$ENV_NAME
  
  if [[ ! -d "$ENV_DIR" ]]; then
      echo "[$(date)] Creating Conda env from $YAML_FILE"
      conda env create --prefix "$ENV_DIR" --file "$YAML_FILE"
  else
      echo "[$(date)] Using existing Conda env at $ENV_DIR"
  fi
  
  # activate env - USE SOURCE not CONDA activate in our environmen
  eval "$(conda shell.bash hook)"
  source activate "$ENV_DIR"
  
  echo "[$(date)] Activated Conda env: $(conda info --envs | grep '*' )"
  
  # run high-memory python task/large matrix operation
  echo "[$(date)] Running high-memory Python script"
  python - <<'EOF'
  import numpy as np, matplotlib.pyplot as plt
  
  # example high-memory operation: large matrix
  arr = np.random.random((200000, 1000))  # ~160 GB
  print("Matrix shape:", arr.shape)
  # Simple plot to use matplotlib
  plt.hist(arr.ravel()[:1000000], bins=50)
  plt.title("Sample histogram")
  plt.savefig("histogram.png")
  print("Histogram saved")
  EOF
  
  echo "[$(date)] Job completed successfully."
  ```
  ```bash
  $ sbatch himem-job.sbatch
  Submitted batch job 238504
  ```
  ```bash
  $ cat matplotlib_238504.out
  [Fri Jun 27 05:54:58 PM MDT 2025] Starting job on host. SLURM_JOB_ID=238504
  [Fri Jun 27 05:55:02 PM MDT 2025] Using existing Conda env at /path/to/conda_envs/matplotlib_3.10.0
  [Fri Jun 27 05:56:00 PM MDT 2025] Activated Conda env:                      * /path/to/conda_envs/matplotlib_3.10.0
  [Fri Jun 27 05:56:01 PM MDT 2025] Running high-memory Python script
  Matrix shape: (200000, 1000)
  Histogram saved
  [Fri Jun 27 05:56:37 PM MDT 2025] Job completed successfully.
  ```

**Explanation**: *Here we use `--partition=himem` (intended for high-memory nodes) and `--mem=256G` to request 256 GB of RAM. You could also use `--mem-per-cpu` to distribute memory per core. Jobs that exceed the memory will fail, so it’s important to request enough. The cluster’s “himem” nodes are designed for these cases.*

---

#### Environment Modules (Lmod)
The cluster uses [Lmod](https://lmod.readthedocs.io/en/latest/) (Lua-based modules) for managing software. Modules let you “load” software packages by adjusting your environment (`PATH`, `LD_LIBRARY_PATH`, etc.) dynamically. Common commands are module avail to list available tools, and `module load <name>` to activate one. For example:
  ```bash
  module avail           # lists all modules - can use “ml” instead of “module”
  module spider samtools
  module load samtools
  which samtools         # points to the samtools executable ( i.e., /data/apps/software/spack/linux-rocky9-x86_64_v3/gcc-11.3.1/samtools-1.21-57fumi7zpuewerxzvu2r3bwkgsuat2gj/bin/samtools )
  ```

In SBATCH scripts, simply include module load lines before running your commands. For instance:
  ```bash
  #!/bin/bash
  #SBATCH --job-name=align
  #SBATCH --output=align_%j.out
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=4
  #SBATCH --mem=8G
  
  module load bwa
  module load samtools
  
  bwa mem ref.fa reads.fq > aln.sam
  ```

Here we load [BWA](https://bio-bwa.sourceforge.net/) and Samtools modules to run an alignment. The Lmod system ensures the correct software versions and dependencies are in your environment.

Apptainer (Singularity) Container Job (with `SIGINT`)
For reproducible or complex environments, [Apptainer](https://apptainer.org/) (formerly Singularity) containers can be used. Apptainer is installed system-wide and used to run containers within your job. For example, an SBATCH script might look like:

  ```bash
  #!/bin/bash
  ##SBATCH --job-name=cleveland_bound
  #SBATCH --output=ctr_%j.out
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=2
  #SBATCH --mem=8G
  
  set -Eeuo pipefail
  trap 'echo "Caught SIGINT in script"; exit 1' INT
  
  srun apptainer exec /path/to/containers/sl.sif sl -e &
  PID=$!
  
  echo "Launched work: PID=$PID; PGID=-$PID"
  
  # Wait a bit, then send simulated Ctrl+C
  sleep 4
  kill -INT -"$PID"
  wait "$PID"
  
        $ sbatch container-job.sbatch
  Submitted batch job 238503
  
  $ cat ctr_238503.out
  
                (@@) (  ) (@)  ( )  @@    ()    @     O     @     O      @
           (   )
       (@@@@)
    (    )
  
  (@@@)
  =        ________                ___________
  __ ___|      _________________
  |  |   H  |  |     |   |         ||_| |_||     _|                \_____A
  |  |   H  |__--------------------| [___] |   =|                        |
  ____|___H__/__|_____/[][]~\_______|       |   -|                        |
  |-----------I_____I [][] []  D   |=======|____|________________________|_
  |=-~O=====O=====O=====O\ ____Y___________|__|__________________________|_
  _|=    ||    ||    ||    |_____/~\___/          |_D__D__D_|  |_D__D__D_|
     \__/  \__/  \__/  \__/      \_/               \_/   \_/    \_/   \_/
  ```

**Explanation**: *This uses `apptainer exec` to run `sl -e` inside of the container `sl.sif` on the allocated resources. This script also replicates `Ctrl+[C]` behavior by capturing the PID of the running apptainer command and using `SIGINT` ( `kill` ) - which can be used to checkpoint/log various points within jobs. You can also use apptainer run if the image has a default runscript. If using GPUs in a container, include the `--nv` flag (e.g., `apptainer exec --nv image.sif` ...) to enable NVIDIA libraries inside the container. Apptainer natively supports Singularity images and commands.*

---

#### Job Arrays, Dependencies, and Workflows
- Job arrays: When you need to run the same task many times (e.g., for different input files or parameters), SLURM job arrays let you bundle them into one submission. Add #SBATCH --array=1-10 to run 10 tasks indexed 1 through 10. Inside the job you can use the environment variable SLURM_ARRAY_TASK_ID to know which array index is running (and pick the corresponding input). You can limit concurrent tasks (e.g. --array=1-100%10 runs at most 10 at a time). Job arrays dramatically simplify batch processing of multiple samples.
  
- **Dependencies**: To chain jobs, use the --dependency flag. For example, one script can sbatch --dependency=afterok:<jobid> next_step.sh so that next_step.sh runs only if <jobid> finished successfully. In practice, users often write a driver script that submits jobs and captures their IDs, then submits dependent jobs. This ensures workflows execute in order (e.g. preprocess → analysis → cleanup).

- **Workflow engines**: Many bioinformatics pipelines use workflow managers like Snakemake or Nextflow, which integrate with SLURM. For instance, Snakemake can submit rules as SLURM jobs by using the --slurm flag (it will handle generating SBATCH calls for each job). You can also write a Nextflow config with process.executor = 'slurm' so that each process is submitted via SBATCH. These tools let you focus on the pipeline logic while leveraging SLURM for parallel execution.

- **Logging and checkpoints**: It’s good practice to write meaningful output and error logs (as above) and periodically save intermediate results in long workflows. That way, if a job fails, you only need to restart from the last successful step rather than from scratch.

- **Profiling**: For performance analysis, you can use Linux tools (time, /usr/bin/time, perf, etc.) inside your job, or SLURM’s accounting (sacct) to see CPU and memory usage. Capturing CPU/GPU utilization over time helps identify bottlenecks. SLURM also supports job steps accounting (sacct -j JOBID --format=...) and sstat (for running jobs) to monitor tasks in real-time.

- **ThinLinc**: For interactive GUI work, the cluster provides a ThinLinc remote desktop service. ThinLinc is a remote X11 desktop solution that offers better performance than standard X forwarding. You connect using the ThinLinc client (over SSH) to get a Linux desktop on the cluster. This allows you to run graphical programs (e.g., RStudio, JupyterLab, IGV, Matlab GUI) as if locally, but executing on HPC nodes. A key advantage is session persistence: you can disconnect and later reconnect without losing open applications.
