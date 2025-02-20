# Using SLURM (Scheduler)

HPC resources are managed by SchedMD's [SLURM workload manager](https://slurm.schedmd.com/quickstart.html). Jobs that will use a significant portion of resources for more than a couple minutes should be launched through SLURM instead of directly on the login/submit nodes. 

!!! info "Video Tutorial"
    An excellent video tutorial to Slurm is available from the Brigham Young University Supercomputing Group. 

    <https://www.youtube.com/watch?v=U42qlYkzP9k>

## Command overview

| Commands |	Syntax |	Description |
| -------- | ------ | ------------ |
| sbatch | sbatch &lt;script_name&gt; | Submit a batch script to Slurm for processing. |
| salloc  | salloc | Opens an interactive shell on the allocated resources |
| squeue | squeue -u | Show information about your job(s) in the queue. The command when run without the -u flag, shows a list of your job(s) and all other jobs in the queue. |
| srun | srun &lt;resource-parameters&gt; | Run jobs interactively on the cluster. |
| skill/scancel | scancel &lt;job-id&gt; | End or cancel a queued job. |
| sacct | sacct | Show information about current and previous jobs. |
| sinfo | sinfo | Get information about the resources on available nodes that make up Skyline. |
| seff | seff &lt;job-id&gt; | Show efficiency of a completed job, e.g. how much memory it actually used vs what was requested, or whether it used all requested CPU cores. |

## Batch Jobs
Batch jobs are the primary way most users will interact with Skyline. A batch job needs a submit script. A script can be as simple as:

``` bash title="test.sbatch"
#!/bin/bash
echo "Hello from $(hostname)"
```
This script is submitted with `sbatch` 

``` bash
sbatch test.sbatch
```

A submit script requires 4 things:

1. **Shebang** The Shebang command tells the shell (which interprets the UNIX commands) to interpret and run the Slurm script using the bash  shell.
    ``` bash
    #!/bin/bash
    ```

1. **Parameters** This section informs Slurm about the name of the job, output filename, amount of RAM, CPU count, nodes, tasks, time, and other parameters to be used for processing the job.

    These SBATCH commands are also know as SBATCH directives and must be preceded with a pound sign and should be in an uppercase format as shown below.

    ``` bash
    #SBATCH --job-name=TestJob
    #SBATCH --output=TestJob.out
    #SBATCH --time=01:00
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem-per-cpu=500M
    ```

1. Load all the software that the project depends on to execute. For example, if you are working on R, you’d definitely require the R software or module to interpret and run your code. Please visit the link [Software](software.md) page for more details about using modules on Skyline.
``` bash
module load r r-bio3d
```
1. Specify the list of tasks to be carried out.
    ```
    Rscript inputfile.R
    ```

``` bash title="submit_script.sbatch"
#!/bin/bash
#SBATCH --job-name=TestJob
#SBATCH --output=TestJob.out
#SBATCH --time=01:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
module load r r-bio3d
Rscript inputfile.R
```

## Launching tasks within jobs
`srun` can be used to launch processes within a job. These tasks are referred to as "job steps" and can be thought of as a sub-allocation; `sbatch` reserves resources and `srun` starts tasks within the job using a subset of those resources. By default, `srun` launches the command in parallel on each CPU core within the job; however, it can also accept many of the same resource targeting arguments as `sbatch`.

!!! Info "Parallel runs"
    srun can also be used as a substitute for mpirun and mpiexec if you load the openmpi module first

## Interactive jobs
If you need to use a compute node interactively (e.g. to debug a GPU program), you can use `salloc` similarly to how you would use `sbatch`, but without a script name--for example, to get an interactive shell on a compute node with one GPU, four CPUs, and 8 GiB memory on one node for 2 hours, you could use: 

``` bash
salloc --time 2:00:00 --nodes 1 --ntasks 4 -p gpu --gres=gpu:1 --mem 8g
```

Once given an allocation, use `srun` to run programs within the allocation

``` bash
srun nvidia-smi
```

`salloc` is handy to run a quick test job before submitting your longer-running job with `sbatch`.

### X11 forwarding for graphical applications

It is possible to run software with a Graphical User Interface (GUI) on Skyline using X11 forwarding. Adding the `-x11` parameter to either `salloc` or `srun` will allow your graphical app to run on the cluster with the same resource requests you would make for any other job, but use a window on your [Thinlinc](https://skyline.niaid.nih.gov/access/#graphical-user-interface) session, as if it's running locally there.

As a simple example, you could use `salloc` to request a single node with 4 CPU cores and 2 GB of RAM and x11 forwarding:

```
salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 -p all --x11 --mem 2g 
```

Or with a GPU:
```
salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --x11 --mem 2g
```

You can then ssh to the node you're assigned with the `-XC` flag, which will establish the X connection, allowing you to then launch your graphical app:

```
$ salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 -pgpu --gres=gpu:1 --x11 --mem 2g
salloc: Granted job allocation 1626951
salloc: Nodes ai-hpcgpu21 are ready for job
$ ssh -XC <node>
$ ml imod/4.11
$ 3dmod
```


## Requesting resources
Resources are requested via parameters using either command line options or in the submit script. [A full list of options is available at on the Slurm website.](https://slurm.schedmd.com/sbatch.html)

Skyline has 2 special resources available. 

- High Memory Nodes

    If you need more than 1TB of memory for any compute job, you can indicate that your job depends on these nodes with the partition resource request modifier `-p himem` or `--partition=himem`
    
    ``` bash
    #!/bin/bash
    #SBATCH --partition=himem
    ./my_program_that_uses_a_lot_of_ram
    ```

- GPUs

    Any job that runs software that relies on a GPU will need to request both the 'gpu' partition (`--partition=gpu`) as well as a number (and optionally type) of GPU (`--gres=gpu:<type>:<number>`)

    ``` bash
    #!/bin/bash
    #SBATCH --partition=gpu
    #SBATCH --gres=gpu:A100:2
    module load namd
    namd test.namd
    ```

## Parameter List
| Directives |   	Description |
| ------- | -------------|
| --job-name=x | Specifies a name for the job allocation. The specified name will appear along with the job id number when querying running jobs on the system. The default is the name of the batch script, or just sbatch if the script is read on sbatch’s standard input. |
| --output=x | Instructs Slurm to connect the batch script’s standard output directly to the filename. If not specified, the default filename is slurm-jobID.out. |
| --partition=x | Requests a specific partition for the resource allocation (gpu, himem, all). If not specified, the default partition is all. |
| --ntasks=x | This option advises the Slurm controller that job steps run within the allocation will launch a maximum of number tasks and offer enough resources. The default is 1 task per node, but note that the --cpus-per-task option will change this default. |
| --cpus-per-task=x | Advises the Slurm controller that ensuing job steps will require ncpus number of processors per task. Without this option, the controller will just try to assign one processor per task. For instance, consider an application that has 4 tasks, each requiring 3 processors. |
| --mem-per-cpu=x | This is the minimum memory required per allocated CPU. Note: It’s highly recommended to specify --mem-per-cpu. If not, the default setting of 500MB will be assigned per CPU. |
| --mail-user=user@mail.net | Defines user who will receive email notification of state changes as defined by --mail-type. |
| --mail-type=x | Notifies user by email when certain event types occur. Valid type values are BEGIN, END, FAIL. The user to be notified is indicated with --mail-user. The values of the --mail-type directive can be declared in one line like so: --mail-type BEGIN, END, FAIL |
| --account=x | Specify which account to run under. This is only needed if you are part of multiple labs or part of a special project to divide time allocation |
| --gres=[type:]x | Requests x GPUs per node, and by default x CPUs.Valid options for [type]: A100 |

<!-- lifted from Kirans documentation in the oeb-ste wiki. -->
## Troubleshooting

### Job Submission Errors
- Invalid account/partition
  
    These are caused due to the fact that the user's does not have a SLURM association in Slurm accounting. Please create a ticket with NIAID OEP HPC Support to request access to the correct partition/association in Slurm.
```
sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
srun: error: Unable to allocate resources: Invalid account or account/partition combination specified
```

### Job Never Starts
- Requested resources not available

    Jobs requesting too many cores or a large amount of memory may not start running very quickly. Run `squeue -u <userid> -t PD` (substitute <userid> with your cluster userID) to see the REASON why jobs are not running. If the REASON seen in `squeue` is `Resources`, then the resources requested in the job submission are not yet available. If the reason is `ReqNodeNotAvail`,  it means a particular node is not yet available either because other jobs are running there, or the node is offline.

- Upcoming maintenance window / reservation

    Prior to a planned outage, the NIAID HPC team will create a [reservation](https://slurm.schedmd.com/reservations.html) to reserve all the cluster nodes during that time, so maintenance can occur. If a job that has not yet started and overlaps the maintenance reservation's start time, it will remain in pending or `PD` state; the job will dispatch after the maintenance period is over. Use `squeue -u <userid>`, will show the  `ReqNodeNotAvail, UnavailableNodes:` as the reason for not starting the jobs.
 
- Dependency has not been satisfied

    The dependency for a job with the `--dependency` parameter that has not been met yet. For example, a job will not start until another job completes. In `squeue -u <userid>`, the REASON why this job is pending is reported as `Dependency`.

- Dependency will never be satisfied

    The dependency for a job with the `--dependency` parameter that will never be met. For example, the job will not start until another job successfully completes. If the previous job failed, the current job will never run. In `squeue`, the REASON why this job is pending is reported as `DependencyNeverSatisfied`. In our current configuration, these jobs will need to be manually removed using `scancel <jobid>`.
 
- Priority

    Jobs with lower priority relative to other jobs in the queue will stay in pending state. In `squeue`, such jobs will have `Priority` reported as the REASON. You can also run `sprio` to see the factors that make up a job's scheduling priority; by default `sprio` will show all jobs in the queue, you can limit this by running sprio -u <userid> (substitute <userid> with your cluster userID).

  ---

## Summary
This document addresses SLURM resource request guidelines for NIAID HPC clusters.

---

## Contents
* [Summary](#Summary)
* [Contents](#Contents)
* [Prerequisites](#Prerequisites)
  * [How many (and which) resources should I request?](#How-many-(and-which)-resources-should-I-request?)
  * [What are the consequences of asking for the wrong resources?](#What-are-the-consequences-of-asking-for-the-wrong-resources?)
* [Procedure](#Procedure)
  * [How do I decide what to request?](#How-do-I-decide-what-to-request?)
    * [Tasks](#The-number-of-tasks)
    * [Nodes](#The-number-of-nodes)
    * [Cores](#The-number-of-cores)
    * [Memory](#The-quantity-of-memory)
    * [Partition](#The-partition)
    * [Node Features with SLURM Constraints](#Requesting-multiple-features/constraints)
  * [Array Jobs](#For-array-jobs)
  * [After the job finishes](#After-the-job-finishes)
* [Conclusion](#Conclusion)


---

## Prerequisites

### How many (and which) resources should I request?

When you submit a job to the SLURM job scheduler on the cluster, you must specify how many resources (nodes, cpus, memory) you need and which partition to submit to. But how to decide? This document will give new users some brief guidance.

*All work done on the cluster needs to be done on resources allocated by SLURM.*

### What are the consequences of asking for the wrong resources?

Why should you care? Why not just copy the same SLURM header from script to script and only worry about it when a job fails?

HPC clusters at NIAID are a community resource. X users submit Y jobs per month. There can sometimes be a significant wait for jobs to begin running. Submitting jobs with inaccurate resource requests can tie things up unnecessarily and cause problems with job scheduling, wait times, and time until completion.

* When you run a job and request too many resources:
  * When the cluster is busy, wait times are unnecessarily extended for everyone. If you only need 1 CPU but request 30, 29 CPUs will sit idle while your job runs.
  * You personally will wait longer for your job to run. To fairly allocate resources among users, job priority is lowered for jobs from users consuming many resources.
* When you run a job and request too few resources:
  * Your job may grind nearly to a halt and take a very long time to complete.
  * Your job may fail, requiring you to resubmit it.
  * The job scheduler may oversubscribe the node you are on, leading to the above problems.

---

## Procedure

### How do I decide what to request?

Ok, now that you're 100% convinced that for the benefit of yourself and the community you need to try your best to accurately tailor your resource requests for each job, how are you supposed to do that?

Well, the answer to this can get complicated, but for most common use-cases in bioinformatics, it really isn't.

To start with, here is a typical slurm header which specifies the details of a job submission:

```
#!/bin/bash 
#SBATCH --job-name=<MY_JOB_NAME>
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=int
#SBATCH --mail-user=<FIRST.LAST>@nih.gov
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
```

There are five lines in there that control the resources requested. You should consider their settings for every job you submit.

```
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=int
```

### The number of tasks

`-n`: The number of tasks you will run. A single number. These are tasks to be scheduled by SLURM and you invoke them using the SLURM command srun. They aren't lines in your script, or anything else like that. This option only applies if you use srun more than once within the script you submit using sbatch. If you don't use srun more than once, or at all, leave this at 1, or leave it out (*the default is 1*).

### The number of nodes

`-N`: The number of nodes you are requesting. You can think of nodes as individual computers that all work together to make up the cluster, and the CPUs as processing units that can operate independently, or be easily harnessed together for harder tasks. For most cases you will want to specify this as **1**.

*To explain*: There are two ways that software packages are commonly parallelized. In the simplest, the software can use many CPUs within a single computer. This is how most 'multithreaded' software is written. If you are running something that has an option for how many 'threads' to use, it refers to this kind of parallelization. In the second, *parallelization occurs on such a large scale that many nodes must coordinate their computation*. This is not common in bioinformatics, and if you are a beginning user, you are unlikely to use software that does this.

So, if you wish to use typical parallelized (*a.k.a. multithreaded*) software, you most likely want to set this number to **1**. If you don't specify it, then the number of CPUs you request (*see below*) may be allocated across nodes. As an example, if you wanted to run an application using 12 cores (a.k.a. threads a.k.a. CPUs) to speed it up, your command might look like this:

```
srun <APP> mem -t 12 reference.fa sample.1.fq sample.2.fq -o sample.sam
```

The `-t 12` indicates you want to use 12 cores. If you don't specify `-N 1`, then your job may be allocated 6 cores on node1 and 6 cores on node2. Then the application will only run on one of those nodes and only use 6 cores. The other 6 cores will sit idle because bwa is not written so that it can coordinate parallel computations across nodes.

### The number of cores

`-c`: The number of cores you are requesting (for our purposes cores and threads are the same thing). This is the number of parallel processes that are going to be run during your job. This number is likely to change for every job you submit. Some programs can use multiple cores simultaneously, and sometimes you may want to run a program multiple times in parallel (using GNU parallel, for example). An example of a program that can use multiple cores is the aligner boa. If you wanted to run bwa with 12 cores (or threads), you would specify `-c 12`. If you fail to specify `-c 12`, then despite telling bwa to use 12 cores, it won't be able to. If you specify `-c 12` in your SLURM header but fail to tell bwa to use 12 cores by using the flag `-t 12` in the bwa command line then bwa will only use 1 core and 11 will sit idle.

**It is critically important that when using multithreaded software, you tell both SLURM and the software how many cores to use.**

How many cores should you request? If you're running a single process, set `-c 12`. *A general rule of thumb is that you should allocate the maximum number that your script will use at once*.

### The quantity of memory

`--mem`: The quantity of memory requested. How much memory does your job need? This is a very tricky question, and it can depend on the particulars of not only the software you are using, but also your data.

Hopefully, the software you are using can give you a sense of how much memory it will need. For example, a genome assembler may tell you that if you have X reads to assemble, then it will require Y memory. This is not always the case, however, and it may not be entirely predictable.

In cases of high uncertainty, the main approach is to make your best guess and then start the job running. After it gets going, you can enter the compute node and use the command top to see how much memory (and CPU) your processes are using. You can do this by first running `squeue`, which will list all your currently running jobs, and the nodes they are running on. Then `srun -w <NODE_NAME> --pty bash` - and finally `top -u <USERNAME>`. 

There aren't hard limits on memory usage, so your process may start using more memory than you asked for. If it is dramatically more, you should cancel the job using `scancel` and resubmit it requesting a higher allocation.

### The partition

`--partition`: The partition you are requesting. You can see the characteristics of the partitions using `sinfo` which lists the sets and their membership. You will notice these sets are overlapping. Most jobs will run on the "int" partition. Some jobs may have specific requirements that require you to use other partitions. One example of this is variant calling with GATK; it is much slower on amd processors. So to avoid having a GATK run assigned to one of these processors you would use a constraint, like `-C ‘rome|intel|skylake’`.

*As defined in the SLURM documentation, when a job requires specific hardware features on a compute node, a user may specfiy feature contstraints in their sbatch file using the `-C` or `--contraint` options*

```
#SBATCH --constraint=<feature list>
```

### Requesting multiple features/constraints

If required, more sopisicated logical constraints may be created using parenthesis and brackets to create grouping of features where parenthesis are used with the AND operator and brackets are use with the OR operator.
Example requesting a amd node with an a100 GPU

```
#!/bin/bash
#SBATCH --job-name=<MY_JOB_NAME>
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=int
#SBATCH --constraint="amd&a100"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<FIRST.LAST>@nih.gov
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
```

### For array jobs

When submitting an array job, ***note that the resources should reflect the needs for a single element of the array, not the entire array***. So if the array was specified as:

```
#SBATCH --array=[1-1000]%20
```

***Indicating 1000 tasks running 20 at a time***, and each task needed 1 processor and 1G memory, you would use `-c 1` and `--mem=1G`, NOT `-c 20` and `--mem=20G`.

---

## Conclusion

### After the job finishes

*Did you request the right resources?*

If you don't use the e-mail option in the SLURM header, you can get this exact same report using the command `seff YOUR_JOBID`.

*e.g.*:

An example of a job requiring relatively few resources would be mapping short reads to a reference genome using bwa mem with 4 cores and 5G memory. Code might look like this,

```
bwa mem -t 4 reference.fa sample.1.fq sample.2.fq -o sample.sam
```

And the requested resources would look like this:

```
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=5G
#SBATCH --partition=int
```

Short read mapping is a low-intensity task because for the most part, each read (or read pair) is handled one at a time.

An example of a job that might require slightly more resources is joint variant calling on a multiple individuals. In particular, memory demands can be highly variable among jobs or genomic regions. The amount of memory required scales with

1. genetic diversity (mostly for standard haplotype and assembly based callers like freebayes and GATK) and

2. the sum of read depth across samples. It scales with genetic diversity because the variant caller needs to simultaneously consider many possible haplotypes, and with read depth because all of the reads overlapping a given position must be held in memory and manipulated simultaneously.

The application Trinity creates de novo transcriptomes from RNA-seq data. Trinity has high memory requirements. It can also use many cores at once (16 is a good number). A typical header would have these options set:

```
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=150G
#SBATCH --partition=int
```

---

<!-- Lifted directly from Jeremy Bell's document on the BigSky Wiki with only very minor edits for Skyline. -->
<!-- https://github.niaid.nih.gov/rmllinux/bigsky/blob/master/documents/slurm-job-arrays.md -->
# Slurm Job Arrays
## Summary
This document outlines the process of submitting job arrays to the SLURM scheduler.

If you want to do something very similar many times, you can split the job up into an array of smaller tasks to parallelize it.

Some examples of situations where array jobs can be helpful:

  * Aligning many fastq files to a reference genome.
  * Evaluating many positions in the genome (e.g. variant calling).
  * Simulating and evaluating many datasets.

---

## Contents
* [Summary](#Summary)
* [Contents](#Contents)
* [Prerequisites](#Prerequisites)
* [Procedure](#Procedure)
  * [Iterating over a list of files](#Iterating-over-a-list-of-files)
  * [Iterating over varying input parameters](#Iterating-over-varying-input-parameters)
* [Conclusion](#Conclusion)
* [See Also](#See-Also)


---

## Prerequisites
An example of a trivial array job script, which would be submitted using `sbatch` is below. To run this as a test, create a directory `mkdir array_test`, save the script there using a text editor, and submit it to the job scheduler using the command `sbatch`.

```shell
#!/bin/bash

#SBATCH --job-name=<JOB_ARRAY_NAME>
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=int
#SBATCH --array=[1-1000]%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<first.last>@nih.gov
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID
```

The job will produce 2 files for each task. Each `.out` file contains the hostname and the task number for each task, and the `.err` file contains any corresponding text printed to the standard error stream.

For each task, the variable `SLURM_ARRAY_TASK_ID` is set to the task number. These numbers can range from 0 to 1000, and can be any set of numbers or ranges, separated by commas. The numbers are specified on this line:

```shell
#SBATCH --array=[1-1000]%20
```

The bracketed range gives the task numbers, and `%20` indicates that 20 tasks should run simultaneously. The other SLURM parameters should be set to match the requirements of each task. For example, if you are aligning a fastq file using `bwa` and you want to use 4 CPUs, you would set `#SBATCH -c 4`.

The key to making these job arrays useful is in how you use the `SLURM_ARRAY_TASK_ID` variable.

---

## Procedure

### Iterating over a list of files
As an example, let's pretend we have Illumina sequencing data from 25 samples. For each sample we have paired sequences in separate files, each named `Sample_A.R1.fastq`, `Sample_A.R2.fastq`, etc. We can make dummy files for this exercise:

```shell
# make a new directory to hold test files and cd into it
mkdir array_test_2
cd array_test_2

touch Sample_{A..Y}.R1.fastq
touch Sample_{A..Y}.R2.fastq
```

Say we wanted to use an array job to align each of these samples in parallel. One way to approach this is to create an array variable containing the list of files to analyze, and then use the `SLURM_ARRAY_TASK_ID` to retrieve elements of the list.

An array job script that does this might look like this:

```shell
#!/bin/bash

#SBATCH --job-name=<JOB_ARRAY_NAME>
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=int
#SBATCH --array=[0-24]%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<first.last>@nih.gov
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *.R1.fastq))

# get specific file name, assign it to FQ1
	# note that FILES variable is 0-indexed so
	# for convenience we also began the task IDs with 0
FQ1=${FILES[$SLURM_ARRAY_TASK_ID]}
# edit the file name to refer to the mate pair file and assign that name to FQ2
FQ2=$(echo $FQ1 | sed 's/R1/R2/')
# create an output file name
OUT=$(echo $FQ1 | sed 's/.R1.fastq/.sam/')
# write the input filenames to the standard output to check that everything ran according to expectations. 
echo $FQ1 $FQ2

echo Files $FQ1 and $FQ2 were aligned by task number $SLURM_ARRAY_TASK_ID on $(date)
```

We created the array variable by `FILES=($(ls -1 *.R1.fastq))`, retrieving all files ending in `*.R1.fastq`.

We then retrieve each individual `FQ1` by using `SLURM_ARRAY_TASK_ID` to grab one element of `FILES`. We then get the paired ends by modifying `FQ1` using `sed: FQ2=$(echo $FQ1 | sed 's/R1/R2/')` and defining the output file, `OUT` similarly.

Another approach to make this slightly more robust would be to generate the FILES array using find:

`FILES=($(find /path/to/input_files/ -name "*R1.fastq"))`

This would use the full path for each file and allow the script to be run from any directory.

If you had a small number of files, you could also avoid using the search pattern altogether and define the array variable inside the script by writing out the file names: `FILES=(Sample_A.R1.fastq, Sample_B.R1.fastq. Sample_C.R1.fastq)`

The approach outlined here would work for any situation where you need to iterate an analysis over many files.

### Iterating over varying input parameters

This approach of defining the list of items to be operated on inside the array job script works fine for lists of files (as long as you define matching patterns, file names, and directories that won't collide). But when the variables you want to iterate over are not straightforward manipulations of the task ID or file names, such as targeting many genomic regions, you need another approach.

In those cases you may want to use a file that contains all the relevant information, and extract what you need for each task.

As an example, let's say we want to call variants on human whole genome sequencing data. Instead of letting one job churn through the whole thing sequentially, we might break the job into 10 megabase chunks and run 20 jobs simultaneously to speed it up.

We can first define these 10mb windows using `bedtools`.

```shell
# make a new directory to hold test files and cd into it. 
mkdir array_test_3
cd array_test_3

# load bedtools
module load bedtools2

# this is a tab delimited file giving chromosome names and lengths for the human genome
	# you can make your own version of this file with a reference genome and "samtools faidx"
GEN=/isg/shared/databases/alignerIndex/animal/hg38_ucsc/hg38_STAR/chrNameLength.txt

bedtools makewindows -g $GEN -w 10000000 >10mb.win.bed
```

The resulting bed file, 10mb.win.bed has 791 lines and looks like this:

```shell
chr1	0	10000000
chr1	10000000	20000000
chr1	20000000	30000000
chr1	30000000	40000000
chr1	40000000	50000000
chr1	50000000	60000000
```

We want to run the variant caller over each region as a separate task. We could specify the job array this way:

```shell
#!/bin/bash

#SBATCH --job-name=<JOB_ARRAY_NAME>
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=int
#SBATCH --array=[1-791]%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<first.last>@nih.gov
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# put relevant data in variables
CHR=$(sed -n ${SLURM_ARRAY_TASK_ID}p 10mb.win.bed | cut -f 1)
START=$(expr $(sed -n ${SLURM_ARRAY_TASK_ID}p 10mb.win.bed | cut -f 2) + 1)
STOP=$(sed -n ${SLURM_ARRAY_TASK_ID}p 10mb.win.bed | cut -f 3)

# define region variable
	# for most tools in the format chr:start-stop
REGION=${CHR}:${START}-${STOP}

echo This task will analyze region $REGION of the human genome. 

# subsequent lines would define input and output file names and the variant caller, e.g:

# module load freebayes
# freebayes -r $REGION -f ref.fa aln.bam >var.${REGION}.vcf
```

We specify the size of the job array as `#SBATCH --array=[1-791]%20`. Note that it begins at 1 in this case. We use `sed` to print only the line of the file corresponding to the array task number and pipe that to `cut` to pull out the column that corresponds to the bit of data we want, placing that in a shell variable. Note that for `START`, we add 1 to the number in the bed file because a quirk of the BED format is that the start column is 0-indexed while the end column is 1-indexed. The region specification format of sequence:start-stop usually expects the start and stop positions to both be 1-indexed.

This would result in 791 jobs being run 20 at a time and produce as many output vcf files. The final step would probably be to use a package like `bcftools` or `vcflib` to filter and combine these into a single vcf.

This approach could be generalized to any situation where you have a repetitive task you wish to parallelize, but input parameters that vary for each task.

---

## See Also
Slurm Job Arrays: https://slurm.schedmd.com/job_array.html

---

### Best Practices 
When submitting your job, it's important to keep in mind these general guidelines:

  * Request only the resources you need and make use of them efficiently.
  * Profile and estimate the actual run time of your job. The cluster efficiency greatly improves with accurate job runtimes. Remember, long run times also means long wait times.
  * Be aware of limits. If you request more resources than the hardware can offer, the scheduler might not reject the job, and it may be stuck in the partition forever.
  * Be aware of the available memory limit. In general, the available memory per core is `(memory_per_node)/(cores_in_use_on_the_node)`.
  * There may be per-user quotas on the system.
  * Send your job to the right partition.
  * Do not run compute-intensive tasks from the submit or the ThinLinc nodes. Doing so slows the nodes, causing login delays for other users and may run into resource limits set by the administrators which can terminate your session, often without notice.
  * It is preferable to run small versus large jobs where possible, for several reasons. First, compared to large jobs, small jobs are easier to schedule and have to wait less time for resources due to backfilling. Second, if your job is terminated unexpectedly, the amount of lost work is minimized. We encourage users to make use of workflow orchestration like [nextflow](https://nextflow.io) and [snakemake](https://snakemake.readthedocs.io/en/stable/). Nextflow and [nf-core](https://nf-co.re) have 100's of curated open-source bioinformatics pipelines and are well-suited for large data workflows.

### Job Management Commands
|Command|Description|
|--|--|
|sacct -j job_id -l | Display job accounting data from a completed job.|
|sbatch script_file | Submit a job.|
|sinfo -o "%20N %16F %8c %9m %15f %17G" | Display a list of available resources.|
|scancel job_id | Delete a job.|
|scontrol hold job_id | Place a job on hold.|
|scontrol release job_id | Release a job from hold.|
|sinfo | Reports the state of queues and nodes managed by Slurm.|
|sinfo -N or --Node | Display a list of nodes|
|scontrol show nodes | Display a list of nodes with detailed node information|
|squeue | Display the list of all jobs across all queues.|
|sacctmgr show qos format="name%-14,Description%-20,priority,maxwall" | Display a neatly formatted list of partition|
|squeue -j job_id | Check the status of a job.|
|squeue -u user_nameor--user=user_name | Check the status of all jobs submitted by the user. Can use $USER in place of user_name.|
|squeue --format "specs" | Display custom job information. For example, :squeue --format"%.10a %.10u %.10q %.12j %.4t %M"|
|sstat job_id | Display information about the resources utilized by a running job or job step. Can only sstat your jobs.|
|seff job_id | Extracts data from the Slurm database and prints a job's efficiency |

Use `man <command>` to learn more and view additional options to pass for any of the commands listed above.

### Job States
When checking the status of a job, the state of a job is listed. Jobs typically pass through several states during their execution. An explanation of the different states are as follows.

|Abbreviation | State | Description |
| -- | -- | -- |
| CA | CANCELLED | The job was explicitly cancelled by the user or system administrator. The job may or may not have been initiated.|
| CD | COMPLETED | The job is completed with an exit code of zero.|
| CG | COMPLETING | The job is in the process of completing.|
| F | FAILED | Job terminated with non-zero exit code or other failure condition. |
| PD | PENDING | The job is queued, eligible to run, or routed.|
| R | RUNNING | The job is running.|
| S | SUSPENDED | Job has an allocation, but execution has been suspended, and CPUs have been released for other jobs.|
| TO | TIMEOUT | Job terminated upon reaching its time limit. |

### Job Reason Codes
In addition to state codes, jobs that are pending will have a “reason code” to explain why the job is pending. Completed jobs will have a reason describing how the job ended. Some codes you might see include:

|Reason | Meaning|
| -- | --|
| AssocMaxJobsLimit | The job is being held because the user/project has hit the limit on running jobs. |
| Dependency | Job has dependencies that have not been met. |
| JobHeldUser | Job is held at user’s request. |
| JobHeldAdmin | Job is held at system administrator’s request. |
| JobLaunchFailure | Job failed to launch (could due to system problems, invalid program name, etc.) |
| NonZeroExitCode | The job exited with some code other than 0 |
| Priority | Other jobs with higher priority exist for the partition/reservation. |
| Reservation | The job is waiting for its reservation to become available. |
| ReqNodeNotAvail | The requested node it’s currently unavailable (it’s in use, reserved, down, draining, etc.) |

---

